"""
Contains the HaplotypeGenerator class that is used to generate simulated
haplotypes and use a Markov Chain to inject pre-calculated
LD patterns into data
"""

import numpy as np
import numpy.typing as npt

from .population import PopulationData


class HaplotypeGenerator:
    """Generates simulated haplotypes from minor allele frequencies and
    joint variant probabilities.

    This class uses a Markov Chain to generate N haplotypes while preserving
    approximate MAF and local + distant LD from adjacent variants.

    Attributes:
       population_data (PopulationData) class from .population containing LD and MAF data
       n_haps (int): The number of haplotypes to simulate
       tau (float): LD signal smoothing hyperparameter, modulates adherence to expected R
       lam (float): decay window, impacts the decay of LD signal as a function of physical distance
       window (int): window size in number of variants around a locus
    """

    tau: float
    lam: float
    window: int
    hap_matrix: npt.NDArray[np.int8]

    def __init__(self, population_data: PopulationData, n_haps: int, **kwargs) -> None:
        """inits HaplotypeGenerator class and sets the initial MAF-weighted haplotypes"""
        self.population_data = population_data
        self.n_haps = n_haps
        self.lam = 7.5
        self.tau = 0.3
        self.window = 2
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.init_haplotype()

    def init_haplotype(self) -> None:
        """pre-sets simulated haplotype matrix based on MAF alone,
        which may be refined with forward/reverse passes"""
        hap_len: int = self.population_data.mafs.shape[0]
        # initialize 2d array. This is dense, but we will store it with int8 (same size as bool)
        # for 100000 variants, that is 1kb per haplotype, which should scale okay
        self.hap_matrix: npt.NDArray[int8] = np.zeros((self.n_haps, hap_len), dtype=np.int8)
        self.hap_matrix[:] = (
            np.random.rand(self.n_haps, hap_len) < self.population_data.mafs[None, :]
        ).astype(int)

    def _get_context(self, a: int, window: int, direction: str) -> list[int]:
        """Gets the context window to the left, right, or both"""
        match direction:
            case "left":
                base = max(0, a - window)
                return list(range(base, a))
            case "right":
                base = min(len(self.population_data.mafs), a + window + 1)
                return list(range(a + 1, base))
            case "both":
                l_base = max(0, a - window)
                r_base = min(len(self.population_data.mafs), a + window + 1)
                return list(range(l_base, a)) + list(range(a + 1, r_base))
            case "_":
                raise ValueError(f"Unknown direction: {direction}")

    def get_prob_in_window(self, a: int, context: list[int]) -> npt.NDArray[np.float64]:
        """For a given window and a target variant, get the P(1|window)"""
        maf_a = self.population_data.mafs[a].clip(1e-9, 1 - 1e-9)
        if len(context) == 0:
            return np.full(self.n_haps, maf_a)
        cs_selected = self.hap_matrix[:, context]
        maf_cs = self.population_data.mafs[context]
        pos_deltas = np.abs(
            self.population_data.positions[context] - self.population_data.positions[a]
        )
        weights = np.exp(-(pos_deltas / self.lam))
        weights /= weights.sum().clip(1e-9, 1 - 1e-9)
        if len(weights) == 0:
            return maf_a
        pa_cs = self.population_data.prob_matrix[a, context].toarray().ravel()
        pa_given_cs = (pa_cs / maf_cs).clip(1e-9, 1 - 1e-9)
        pa_given_ncs = ((maf_a - pa_cs) / (1.0 - maf_cs)).clip(1e-9, 1 - 1e-9)
        # get the final p for each marker by multiplying selected by cs and not selected by ncs
        ps = np.where(cs_selected == 1, pa_given_cs, pa_given_ncs).clip(1e-9, 1 - 1e-9)
        lp1 = (np.log(ps) * weights).sum(axis=1) * self.tau + np.log(maf_a)
        lp0 = (np.log((1.0 - ps).clip(1e-9, 1 - 1e-9)) * weights).sum(
            axis=1
        ) * self.tau + np.log(1.0 - maf_a)
        p = 1.0 / (1.0 + (np.exp(lp0 - lp1)))
        return p

    def forward_pass(self) -> None:
        """Do a left to right pass, only looking backward (to the left)"""
        for i in range(1, self.hap_matrix.shape[1] - 1):
            context = self._get_context(i, self.window, direction="left")
            p = self.get_prob_in_window(i, context)
            chosen: npt.NDArray[np.int_] = (np.random.rand(self.n_haps) < p).astype(int)
            self.hap_matrix[:, i] = chosen

    def reverse_pass(self) -> None:
        """Do a right to left pass, looking to the left and right of each target"""
        for i in reversed(range(0, self.hap_matrix.shape[1] - 1)):
            context = self._get_context(i, self.window, direction="both")
            p = self.get_prob_in_window(i, context)
            chosen: npt.NDArray[np.int_] = (np.random.rand(self.n_haps) < p).astype(int)
            self.hap_matrix[:, i] = chosen
