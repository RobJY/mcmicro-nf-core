/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import groovy.io.FileType
import nextflow.Nextflow

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mcmicro_pipeline'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BASICPY                } from '../modules/nf-core/basicpy/main'
include { ASHLAR                 } from '../modules/nf-core/ashlar/main'
include { BACKSUB                } from '../modules/nf-core/backsub/main'
include { CELLPOSE               } from '../modules/nf-core/cellpose/main'
include { DEEPCELL_MESMER        } from '../modules/nf-core/deepcell/mesmer/main'
include { MCQUANT                } from '../modules/nf-core/mcquant/main'
include { SCIMAP_MCMICRO         } from '../modules/nf-core/scimap/mcmicro/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MCMICRO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    input_type = params.input_cycle ? "cycle" : "sample"

    // Run Illumination Correction
    if ( params.illumination ) {

        if (params.illumination == 'basicpy') {

            ch_samplesheet
                .transpose()
                .map { [[it[1].split('/')[-1][0..-5],it[0]], it[1]] }
                .set { ashlar_input_keyed }

                ch_samplesheet
                    .transpose()
                    .set { ch_basicpy_input }

            BASICPY(ch_basicpy_input)
            ch_versions = ch_versions.mix(BASICPY.out.versions)

            BASICPY.out.fields
                .transpose()
                .map { [[it[1].getBaseName()[0..-5],it[0]], it[1]]}
                .groupTuple()
                .set { correction_files_keyed }

            ashlar_input_keyed
                .concat(correction_files_keyed)
                .groupTuple()
                .map { [it[0][1], it[1][1]] }
                .transpose()
                .branch {
                    dfp: it =~ /-dfp.tiff/
                    ffp: it =~ /-ffp.tiff/
                }
                .set { ordered_correction_files }
            ch_dfp = ordered_correction_files.dfp
                .groupTuple()
                .map { it[1] }
            ch_ffp = ordered_correction_files.ffp
                .groupTuple()
                .map { it[1] }

        } else if(params.illumination == 'manual') {

            if (input_type == "cycle") {
                samplesheet = "input_cycle"
                dfp_index = 4
                ffp_index = 5
            } else if (input_type == "sample") {
                samplesheet = "input_sample"
                dfp_index = 3
                ffp_index = 4
            }
            ch_manual_illumination_correction = Channel.fromSamplesheet(
                samplesheet,
                skip_duplicate_check: false
            )
            .multiMap
                { it ->
                    dfp: it[dfp_index]
                    ffp: it[ffp_index]
                }

            ch_dfp = ch_manual_illumination_correction.dfp
            ch_ffp = ch_manual_illumination_correction.ffp
        }

    } else {
        ch_dfp = []
        ch_ffp = []
    }

    // Run Registration
    ASHLAR(ch_samplesheet.map{[[id:it[0]['id']], it[1][0]]}, ch_dfp, ch_ffp)
    ch_versions = ch_versions.mix(ASHLAR.out.versions)

    // Run Background Correction
    if (params.backsub) {
        BACKSUB(ASHLAR.out.tif, [[id:"$ASHLAR.out.tif[0]['id']"], params.marker_sheet])
        ch_versions = ch_versions.mix(BACKSUB.out.versions)
    }

    // Run Segmentation
    if (params.segmentation == "mesmer") {
        DEEPCELL_MESMER(ASHLAR.out.tif, [[:],[]])
        ch_versions = ch_versions.mix(DEEPCELL_MESMER.out.versions)
        mcquant_in = ASHLAR.out.tif.join(DEEPCELL_MESMER.out.mask).multiMap { it ->
            image: [it[0], it[1]]
            mask: [it[0], it[2]]
        }
    } else if (params.segmentation = "cellpose") {
        CELLPOSE( ASHLAR.out.tif, [] )
        ch_versions = ch_versions.mix(CELLPOSE.out.versions)
        mcquant_in = ASHLAR.out.tif.join(CELLPOSE.out.mask).multiMap { it ->
            image: [it[0], it[1]]
            mask: [it[0], it[2]]
        }
    } else if (params.segmentation = "unmicst"){
        error("apologies, unmicst not supported yet")
    }

    // Run Quantification
    MCQUANT(mcquant_in.image,
            mcquant_in.mask,
            [[:], file(params.marker_sheet)])
    ch_versions = ch_versions.mix(MCQUANT.out.versions)

    // Run Reporting
    SCIMAP_MCMICRO(MCQUANT.out.csv)
    ch_versions = ch_versions.mix(SCIMAP_MCMICRO.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
