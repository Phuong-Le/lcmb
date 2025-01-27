include { getPhylogeny } from "$projectDir/modules/local/get_phylogeny"
include { matrixGeneratorBranches } from "$projectDir/modules/local/matrix_generator_branches"
include { concatMatrices } from "$projectDir/modules/local/concat_matrices"
include { spectraPlottingBranches } from "$projectDir/modules/local/spectra_plotting_branches"


workflow PHYLOGENETICS { // phylogenetics workflow for SNVs
    take:
    phylogenetics_input_ch
    outdir_basename
    sigprofiler_genome

    main:
    //  get phylogeny
    getPhylogeny(
        phylogenetics_input_ch,
        outdir_basename
    )

    // generate mutation matrix for the branches by SigProfilerMatrixGenerator
    matrixGeneratorBranches(
        getPhylogeny.out.muts_on_tree,
        outdir_basename,
        sigprofiler_genome
        )

    concatMatrices(
        matrixGeneratorBranches.out.toList(),
        outdir_basename
        )

    // // plotting
    // spectraPlottingBranches(
    //     concatMatrices.out,
    //     outdir_basename
    //     )

    emit:
    getPhylogeny.out.topology

}
