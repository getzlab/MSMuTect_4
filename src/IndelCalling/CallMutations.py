# cython: language_level=3
import numpy as np
from scipy.stats import binom

from src.IndelCalling.CallAlleles import calculate_alleles
from src.IndelCalling.FisherTest import Fisher
from src.IndelCalling.MutationCall import MutationCall
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.AICs import AICs

# used for generating sets for fisher test.
# ex. (5.0_4, 6.0_5) and (3.0_2, 5.0_1) -> first_set = [0, 4, 5], second_set = [2, 0, 1]
from src.IndelCalling.hist2vecs import hist2vecs


def cdf_test(first_allele_reads: int, second_allele_reads: int, p_equal: float = 0.3):
    p = binom.cdf(min(first_allele_reads, second_allele_reads), first_allele_reads + second_allele_reads, 0.5)
    if p < p_equal:
        return MutationCall.INSUFFICIENT
    else:
        return MutationCall.MUTATION


def check_normal_alleles(normal_alleles: AlleleSet, p_equal=0.3) -> int:
    if len(normal_alleles.repeat_lengths) == 0:
        return MutationCall.NO_ALLELES
    elif len(normal_alleles.repeat_lengths) == 1:
        return MutationCall.MUTATION  # not that it is necessarily a mutation yet; however it may be
    elif len(normal_alleles.repeat_lengths) == 2:
        first_allele = normal_alleles.repeat_lengths[0]
        second_allele = normal_alleles.repeat_lengths[1]
        first_allele_reads = normal_alleles.histogram.repeat_lengths[first_allele]
        second_allele_reads = normal_alleles.histogram.repeat_lengths[second_allele]
        return cdf_test(first_allele_reads, second_allele_reads, p_equal)
    else: 
        return MutationCall.TOO_MANY_ALLELES


def log_likelihood(histogram: Histogram, alleles: AlleleSet, noise_table: np.array) -> float:
    L_k_log = 0
    rounded_histogram = histogram.rounded_repeat_lengths
    if alleles.repeat_lengths.size == 0:
        return -1_000_000.0
    for length in rounded_histogram.keys():
        if length < 40:
            L_k_log+=rounded_histogram[length]*np.log(sum(alleles.frequencies*noise_table[alleles.repeat_lengths, length]) + 1e-6)
    return L_k_log


def calculate_AICs(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table: np.array) -> AICs:
    num_normal_alleles = len(normal_alleles.repeat_lengths)
    num_tumor_alleles = len(tumor_alleles.repeat_lengths)
    L_Norm_Tum = log_likelihood(normal_alleles.histogram, tumor_alleles, noise_table)
    L_Norm_Norm = log_likelihood(normal_alleles.histogram, normal_alleles, noise_table)
    L_Tum_Tum = log_likelihood(tumor_alleles.histogram, tumor_alleles, noise_table)
    L_Tum_Norm = log_likelihood(tumor_alleles.histogram, normal_alleles, noise_table)
    return AICs(normal_normal=2 * num_normal_alleles - 2 * L_Norm_Norm,
                normal_tumor=2 * num_tumor_alleles - 2 * L_Norm_Tum,
                tumor_tumor=2 * num_tumor_alleles - 2 * L_Tum_Tum,
                tumor_normal=2 * num_normal_alleles - 2 * L_Tum_Norm)


def passes_AICs(AIC_scores: AICs, LOR_ratio = 8.0) -> bool:
    return (AIC_scores.tumor_tumor - AIC_scores.tumor_normal) < -LOR_ratio \
    and (AIC_scores.normal_normal - AIC_scores.normal_tumor) < -LOR_ratio


def fisher_test(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, fisher_calculator: Fisher) -> float:
    reads_sets = hist2vecs(tumor_alleles.histogram, normal_alleles.histogram)  # order is important for Fisher test
    one_sided_fisher = fisher_calculator.test(reads_sets.first_set, reads_sets.second_set)
    return one_sided_fisher


def call_decision(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table: np.array, fisher_calculator,
                  LOR_ratio = 8.0, p_equal = 0.3, fisher_threshold = 0.031) -> MutationCall:
    normal_allele_call = check_normal_alleles(normal_alleles, p_equal)
    if normal_allele_call != MutationCall.MUTATION:
        return MutationCall(normal_allele_call, normal_alleles, tumor_alleles, AICs())
    else:
        return call_verified_locus(normal_alleles, tumor_alleles, noise_table, fisher_calculator, fisher_threshold, LOR_ratio)


def equivalent_arrays(a: np.array, b: np.array) -> bool:
    # returns whether arrays hold the same values in any value, is intended for integer arrays
    if len(a)!=len(b):
        return False
    a_sorted = np.sort(a.astype(np.int32))
    b_sorted = np.sort(b.astype(np.int32))
    return (a_sorted == b_sorted).all()


def call_mutations(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table: np.array, fisher_calculator: Fisher) -> MutationCall:
    if len(normal_alleles) == 0 or len(tumor_alleles) == 0:
        return MutationCall(MutationCall.NO_ALLELES, normal_alleles, tumor_alleles, AICs())
    elif equivalent_arrays(normal_alleles.repeat_lengths, tumor_alleles.repeat_lengths):
        return MutationCall(MutationCall.NOT_MUTATION, normal_alleles, tumor_alleles, AICs())
    else:
        return call_decision(normal_alleles, tumor_alleles, noise_table, fisher_calculator)


def is_possible_mutation(normal_alleles: AlleleSet, p_equal = 0.3) -> bool:
    # checks if locus is a candidate to be called a mutation based on its normal alleles
    return check_normal_alleles(normal_alleles, p_equal) == MutationCall.MUTATION


def reconstruct_tumor_alleles_without_reference_length(tumor_alleles: AlleleSet, noise_table, integer_indels_only: bool) -> AlleleSet:
        new_locus = tumor_alleles.histogram.locus
        ref_length = new_locus.repeats
        new_histogram = Histogram(new_locus, integer_indels_only)
        new_histo_dict = dict()
        for repeat in tumor_alleles.histogram.rounded_repeat_lengths.keys():
            new_histo_dict[repeat] = tumor_alleles.histogram.rounded_repeat_lengths[repeat]
        if ref_length in new_histo_dict.keys():
            del new_histo_dict[ref_length]
        new_histogram.repeat_lengths = new_histo_dict
        new_tumor_alleles = calculate_alleles(new_histogram, noise_table, required_read_support=tumor_alleles.min_read_support)
        return new_tumor_alleles


def reversion_to_reference_simple(tumor_alleles: AlleleSet) -> bool:
    reference_length = int(tumor_alleles.histogram.locus.repeats)
    return len(tumor_alleles.repeat_lengths)==1 and reference_length in tumor_alleles.repeat_lengths.astype(np.int32)


def reversion_to_reference(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table: np.array, fisher_calculator: Fisher,
                           fisher_threshold = 0.031, LOR_ratio = 8.0) -> bool:
    reference_length = normal_alleles.histogram.locus.repeats
    if reference_length not in tumor_alleles.repeat_lengths:
        return False
    tumor_alleles_ref_removed = reconstruct_tumor_alleles_without_reference_length(tumor_alleles, noise_table, tumor_alleles.histogram.integer_indels_only)
    if equivalent_arrays(normal_alleles.repeat_lengths, tumor_alleles_ref_removed.repeat_lengths):
        return True
    aic_values = calculate_AICs(normal_alleles, tumor_alleles_ref_removed, noise_table)
    if passes_AICs(aic_values, LOR_ratio):
        p_value = fisher_test(normal_alleles, tumor_alleles_ref_removed, fisher_calculator)
        if p_value < fisher_threshold:
            return False
        else:
            return True
    else:
        return True


def call_verified_locus(normal_alleles: AlleleSet, tumor_alleles: AlleleSet, noise_table: np.array, fisher_calculator: Fisher,
                        fisher_threshold = 0.031, LOR_ratio = 8.0) -> MutationCall:
    # calls mutation for locus that has proper normal alleles and support
    aic_values = calculate_AICs(normal_alleles, tumor_alleles, noise_table)
    if passes_AICs(aic_values, LOR_ratio):
        p_value = fisher_test(normal_alleles, tumor_alleles, fisher_calculator)
        if p_value < fisher_threshold:
            if reversion_to_reference(normal_alleles, tumor_alleles, noise_table, fisher_calculator, fisher_threshold, LOR_ratio):
                return MutationCall(MutationCall.REVERTED_TO_REFERENCE, normal_alleles, tumor_alleles, aic_values, p_value)
            else:
                return MutationCall(MutationCall.MUTATION, normal_alleles, tumor_alleles, aic_values, p_value)
        else:
            return MutationCall(MutationCall.BORDERLINE_NONMUTATION, normal_alleles, tumor_alleles, aic_values, p_value)
    else:
        return MutationCall(MutationCall.NOT_MUTATION, normal_alleles, tumor_alleles, aic_values)
