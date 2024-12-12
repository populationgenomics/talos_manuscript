import zoneinfo
from datetime import datetime
from enum import Enum

from pydantic import BaseModel, Field

# Steal the talos models

NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
CHROM_ORDER = list(map(str, range(1, 23))) + NON_HOM_CHROM

# some kind of version tracking
CURRENT_VERSION = '1.1.0'
ALL_VERSIONS = [None, '1.0.0', '1.0.1', '1.0.2', '1.0.3', '1.1.0']

# ratios for use in AB testing
MAX_WT = 0.15
MIN_HET = 0.25
MAX_HET = 0.75
MIN_HOM = 0.85
_GRANULAR_DATE = None
TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')


def get_granular_date():
    """
    cached getter/setter
    """
    global _GRANULAR_DATE
    if _GRANULAR_DATE is None:
        _GRANULAR_DATE = datetime.now(tz=TIMEZONE).strftime('%Y-%m-%d')
    return _GRANULAR_DATE

class FileTypes(Enum):
    """
    enumeration of permitted input file types
    """

    HAIL_TABLE = '.ht'
    MATRIX_TABLE = '.mt'
    PED = 'ped'
    VCF = '.vcf'
    VCF_GZ = '.vcf.gz'
    VCF_BGZ = '.vcf.bgz'


class PhenoPacketHpo(BaseModel):
    """
    A representation of a HPO term
    """

    id: str
    label: str


class Coordinates(BaseModel):
    """
    A representation of genomic coordinates
    """

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def string_format(self) -> str:
        """
        forms a string representation: chr-pos-ref-alt
        """
        return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'

    def __lt__(self, other) -> bool:
        """
        enables positional sorting
        """
        # this will return False for same chrom and position
        if self.chrom == other.chrom:
            return self.pos < other.pos
        # otherwise take the relative index from sorted chromosomes list
        if self.chrom in CHROM_ORDER and other.chrom in CHROM_ORDER:
            return CHROM_ORDER.index(self.chrom) < CHROM_ORDER.index(other.chrom)
        # if self is on a canonical chromosome, sort before HLA/Decoy etc.
        if self.chrom in CHROM_ORDER:
            return True
        return False


class VariantCommon(BaseModel):
    """
    the abstracted representation of a variant from any source
    """

    coordinates: Coordinates = Field(repr=True)
    info: dict[str, str | int | float | list[str] | list[float] | dict[str, str] | bool] = Field(default_factory=dict)
    het_samples: set[str] = Field(default_factory=set, exclude=True)
    hom_samples: set[str] = Field(default_factory=set, exclude=True)
    boolean_categories: list[str] = Field(default_factory=list, exclude=True)
    sample_categories: list[str] = Field(default_factory=list, exclude=True)
    sample_support: list[str] = Field(default_factory=list, exclude=True)
    phased: dict = Field(default_factory=dict)

    def __str__(self):
        return repr(self)

    def __lt__(self, other):
        return self.coordinates < other.coordinates

    def __eq__(self, other):
        return self.coordinates == other.coordinates

    @property
    def has_boolean_categories(self) -> bool:
        """
        check that the variant has at least one assigned class
        """
        return any(self.info[value] for value in self.boolean_categories)

    @property
    def has_sample_categories(self) -> bool:
        """
        check that the variant has any list-category entries
        """
        return any(self.info[value] for value in self.sample_categories)

    @property
    def has_support(self) -> bool:
        """
        check for a True flag in any CategorySupport* attribute
        Returns:
            True if variant is support
        """
        return any(self.info[value] for value in self.sample_support)

    @property
    def category_non_support(self) -> bool:
        """
        check the variant has at least one non-support category assigned
        Returns:
            True if var has a non-support category assigned
        """
        return self.has_sample_categories or self.has_boolean_categories

    @property
    def is_classified(self) -> bool:
        """
        check for at least one assigned class, inc. support
        Returns:
            True if classified
        """
        return self.category_non_support or self.has_support

    @property
    def support_only(self) -> bool:
        """
        check that the variant is exclusively cat. support
        Returns:
            True if support only
        """
        return self.has_support and not self.category_non_support

    def sample_support_only(self, sample_id: str) -> bool:
        """
        check that the variant is exclusively cat. support
        check that this sample is missing from sample flags

        Returns:
            True if support only
        """
        return self.has_support and not (self.category_non_support or self.sample_categorised_check(sample_id))

    def category_values(self, sample: str) -> set[str]:
        """
        get all variant categories
        steps category flags down to booleans - true for this sample

        Args:
            sample (str): sample id

        Returns:
            set of all categories applied to this variant
        """

        # step down all category flags to boolean flags
        categories: set[str] = set()
        for category in self.sample_categories:
            cat_samples = self.info[category]
            if not isinstance(cat_samples, list):
                raise TypeError(f'Sample categories should be a list: {cat_samples}')
            if sample in cat_samples:
                categories.add(category.removeprefix('categorysample'))

        categories.update(
            {bool_cat.replace('categoryboolean', '') for bool_cat in self.boolean_categories if self.info[bool_cat]},
        )

        if self.has_support:
            categories.add('support')

        return categories

    def sample_categorised_check(self, sample_id: str) -> bool:
        """
        check if any *sample categories applied for this sample

        Args:
            sample_id (str): the specific sample ID to check

        Returns:
            bool: True if this sample features in any
                  named-sample category, includes 'all'
        """
        for category in self.sample_categories:
            cat_samples = self.info[category]
            assert isinstance(cat_samples, list)
            if any(sam in cat_samples for sam in [sample_id, 'all']):
                return True

        return False

    def sample_category_check(self, sample_id: str, allow_support: bool = True) -> bool:
        """
        take a specific sample and check for assigned categories
        optionally, include checks for support category

        Args:
            sample_id (str):
            allow_support: (bool) also check for support

        Returns:
            True if the variant is categorised for this sample
        """
        big_cat = self.has_boolean_categories or self.sample_categorised_check(sample_id)
        if allow_support:
            return big_cat or self.has_support
        return big_cat

    def check_ab_ratio(self, *args, **kwargs) -> set[str]:  # noqa: ARG002, ANN002, ANN003
        """
        dummy method for AB ratio checking - not implemented for SVs
        """
        return set()

    def get_sample_flags(self, *args, **kwargs) -> set[str]:  # noqa: ARG002, ANN002, ANN003
        """
        dummy method for flag checking - not implemented for SVs (yet)
        """
        return set()

    def check_read_depth(self, *args, **kwargs) -> set[str]:  # noqa: ARG002, ANN002, ANN003
        """
        dummy method for read depth checking - not implemented for SVs
        """
        return set()


class SmallVariant(VariantCommon):
    """
    a representation of a small variant
    note that transcript_consequences is not optional
    we require that something specific to SmallVariant(s) is mandatory
    this is in order to correctly deserialise the Small/Structural objects
    into the appropriate types. If it is optional, Pydantic can coerce
    everything as a SmallVariant
    """

    depths: dict[str, int] = Field(default_factory=dict, exclude=True)
    ab_ratios: dict[str, float] = Field(default_factory=dict, exclude=True)
    transcript_consequences: list[dict[str, str | float | int]]

    def get_sample_flags(self, sample: str) -> set[str]:
        """
        gets all report flags for this sample - currently only one flag
        """
        return self.check_ab_ratio(sample) | self.check_read_depth(sample)

    def check_read_depth(self, sample: str, threshold: int = 10, var_is_cat_1: bool = False) -> set[str]:
        """
        flag low read depth for this sample

        Args:
            sample (str): sample ID matching VCF
            threshold (int): cut-off for flagging as a failure
            var_is_cat_1 (bool): flag if this variant is a category 1

        Returns:
            return a flag if this sample has low read depth
        """
        if self.depths[sample] < threshold and not var_is_cat_1:
            return {'Low Read Depth'}
        return set()

    def check_ab_ratio(self, sample: str) -> set[str]:
        """
        AB ratio test for this sample's variant call

        Args:
            sample (str): sample ID

        Returns:
            set[str]: empty, or indicating an AB ratio failure
        """
        het = sample in self.het_samples
        hom = sample in self.hom_samples
        variant_ab = self.ab_ratios.get(sample, 0.0)

        if (variant_ab <= MAX_WT) or (het and not MIN_HET <= variant_ab <= MAX_HET) or (hom and variant_ab <= MIN_HOM):
            return {'AB Ratio'}
        return set()


class StructuralVariant(VariantCommon):
    """
    placeholder for any methods/data specific to Structural Variants
    """


# register all interchangeable models here
VARIANT_MODELS = SmallVariant | StructuralVariant


class ReportPanel(BaseModel):
    """
    simple storage for all the panels to present in tooltips
    """

    forced: dict[int, str] = Field(default_factory=dict)
    matched: dict[int, str] = Field(default_factory=dict)


class ReportVariant(BaseModel):
    """
    A variant passing MOI tests, to be reported
    """

    sample: str
    var_data: VARIANT_MODELS
    categories: set[str] = Field(default_factory=set)
    date_of_phenotype_match: str | None = None
    phenotype_labels: set[str] = Field(default_factory=set)

    evidence_last_updated: str = Field(default=get_granular_date())
    family: str = Field(default_factory=str)
    # 'tagged' is seqr-compliant language
    first_tagged: str = Field(default=get_granular_date())
    flags: set[str] = Field(default_factory=set)
    gene: str = Field(default_factory=str)
    genotypes: dict[str, str] = Field(default_factory=dict)
    independent: bool = Field(default=False)
    labels: set[str] = Field(default_factory=set)
    panels: ReportPanel = Field(default_factory=ReportPanel)
    phenotypes: list[PhenoPacketHpo] = Field(default_factory=list)
    reasons: set[str] = Field(default_factory=set)
    support_vars: set[str] = Field(default_factory=set)
    # log whether there was an increase in ClinVar star rating since the last run
    clinvar_increase: bool = Field(default=False)

    def __eq__(self, other):
        """
        makes reported variants comparable
        """
        return self.sample == other.sample and self.var_data.coordinates == other.var_data.coordinates

    def __lt__(self, other):
        return self.var_data.coordinates < other.var_data.coordinates


class PanelDetail(BaseModel):
    """
    A gene from PanelApp, combining all MOI and panel IDs
    where the gene features on multiple panels
    """

    symbol: str
    chrom: str = Field(default_factory=str)
    all_moi: set[str] = Field(default_factory=set)
    moi: str = Field(default_factory=str)
    new: set[int] = Field(default_factory=set)
    panels: set[int] = Field(default_factory=set)


class PanelShort(BaseModel):
    """
    Short panel summary, used in the metadata section
    """

    id: int
    name: str = Field(default_factory=str)
    version: str = 'UNKNOWN'


class PanelApp(BaseModel):
    metadata: list[PanelShort] = Field(default_factory=list)
    genes: dict[str, PanelDetail]
    version: str = CURRENT_VERSION


class CategoryMeta(BaseModel):
    """
    The mapping of category names to their display names
    """

    categories: dict[str, str] = Field(default=dict)


class HistoricSampleVariant(BaseModel):
    """ """

    # categories here will be a dict of {categories: associated date first seen}
    categories: dict[str, str]
    # new variable to store the date the variant was first seen, static
    first_tagged: str
    support_vars: set[str] = Field(
        default_factory=set,
        description='supporting variants if this has been identified in a comp-het',
    )
    independent: bool = Field(default=True)
    clinvar_stars: int | None = None
    first_phenotype_tagged: str | None = None
    phenotype_labels: set[str] = Field(default_factory=set)


class HistoricVariants(BaseModel):
    """
    The model representing the state transition file
    All relevant metadata relating to the available categories
    Then a per-participant dict of variants, containing the categories
    they have been assigned, date first seen, and supporting variants
    """

    metadata: CategoryMeta = Field(default_factory=CategoryMeta)
    # dict - participant ID -> variant -> variant data
    results: dict[str, dict[str, HistoricSampleVariant]] = Field(default_factory=dict)
    version: str = CURRENT_VERSION


class ResultMeta(BaseModel):
    """
    metadata for a result set
    """

    categories: dict[str, str] = Field(default=dict)
    version: str = Field(default_factory=str)
    family_breakdown: dict[str, int] = Field(default_factory=dict)
    input_file: str = Field(default_factory=str)
    panels: list[PanelShort] = Field(default_factory=list)
    run_datetime: str = Field(default=get_granular_date())
    projects: list[str] = Field(default_factory=list)


class MemberSex(Enum):
    MALE = 'male'
    FEMALE = 'female'
    UNKNOWN = 'unknown'


class FamilyMembers(BaseModel):
    affected: bool = Field(default=False)
    ext_id: str = Field(default_factory=str)
    sex: MemberSex = Field(default=MemberSex.UNKNOWN)


class ParticipantMeta(BaseModel):
    ext_id: str
    family_id: str
    members: dict[str, FamilyMembers] = Field(default_factory=dict)
    phenotypes: list[PhenoPacketHpo] = Field(default_factory=list)
    panel_details: dict[int, str] = Field(default_factory=dict)
    solved: bool = Field(default=False)
    present_in_small: bool = Field(default=False)
    present_in_sv: bool = Field(default=False)


class ParticipantResults(BaseModel):
    """
    A representation of a result set
    """

    variants: list[ReportVariant] = Field(default_factory=list)
    metadata: ParticipantMeta = Field(default_factory=ParticipantMeta)


class ResultData(BaseModel):
    """
    A representation of a result set
    """

    results: dict[str, ParticipantResults] = Field(default_factory=dict)
    metadata: ResultMeta = Field(default_factory=ResultMeta)
    version: str = CURRENT_VERSION


class ModelVariant(BaseModel):
    """
    might be required for the VCF generator
    """


class ParticipantHPOPanels(BaseModel):
    external_id: str = Field(default_factory=str)
    family_id: str = Field(default_factory=str)
    hpo_terms: list[PhenoPacketHpo] = Field(default_factory=list)
    panels: set[int] = Field(default_factory=set)
    matched_genes: set[str] = Field(default_factory=set)
    matched_phenotypes: set[str] = Field(default_factory=set)


class PhenotypeMatchedPanels(BaseModel):
    samples: dict[str, ParticipantHPOPanels] = Field(default_factory=dict)
    all_panels: set[int] = Field(default_factory=set)
    version: str = CURRENT_VERSION


class MiniVariant(BaseModel):
    categories: set[str] = Field(default_factory=set)
    support_vars: set[str] = Field(default_factory=set)
    independent: bool = Field(default=True)


class MiniForSeqr(BaseModel):
    metadata: CategoryMeta = Field(default_factory=CategoryMeta)
    results: dict[str, dict[str, MiniVariant]] = Field(default_factory=dict)


class PedigreeMember(BaseModel):
    """
    This will be a more searchable implementation of the peds pedigree
    """

    family: str
    id: str
    mother: str | None = None
    father: str | None = None
    sex: str
    affected: str
    ext_id: str = 'Missing'
    hpo_terms: list[PhenoPacketHpo] = Field(default_factory=list)


class Pedigree(BaseModel):
    members: list[PedigreeMember] = Field(default_factory=list)
    by_family: dict[str, list[PedigreeMember]] = Field(default_factory=dict)
    by_id: dict[str, PedigreeMember] = Field(default_factory=dict)
