/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import groovy.io.FileType
import nextflow.Nextflow

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

//WorkflowMcmicro.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { BASICPY                     } from '../modules/nf-core/basicpy/main'
include { ASHLAR                      } from '../modules/nf-core/ashlar/main'
include { BACKSUB                     } from '../modules/nf-core/backsub/main'
include { CELLPOSE                    } from '../modules/nf-core/cellpose/main'
include { DEEPCELL_MESMER             } from '../modules/nf-core/deepcell/mesmer/main'
include { MCQUANT                     } from '../modules/nf-core/mcquant/main'
include { SCIMAP_MCMICRO              } from '../modules/nf-core/scimap/mcmicro/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { ILASTIK_PIXELCLASSIFICATION } from '../modules/nf-core/ilastik/pixelclassification/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Manually define inputs here
//image_tuple = tuple([ id:'image' ], '/home/florian/Documents/tmp_data_folder/cycif_tonsil_registered.ome.tif')
//marker_tuple = tuple([ id:'marker'], '/home/florian/Documents/tmp_data_folder/markers.csv')

// Info required for completion email and summary
def multiqc_report = []

workflow MCMICRO {

    def input_type

    if (params.input_sample && !params.input_cycle) {
        input_type = "sample"
        sample_sheet_index_map = make_sample_sheet_index_map(params.input_sample)
        ch_from_samplesheet = Channel.fromSamplesheet("input_sample")
            .multiMap
                { it ->
                    ashlar: make_ashlar_input_sample(it, sample_sheet_index_map)
                }
    } else if(!params.input_sample && params.input_cycle) {
        input_type = "cycle"
        sample_sheet_index_map = make_sample_sheet_index_map(params.input_cycle)
        ch_from_samplesheet = Channel.fromSamplesheet("input_cycle")
            .multiMap
                { it ->
                    ashlar: make_ashlar_input_cycle(it, sample_sheet_index_map)
                }
    } else if(params.input_sample && params.input_cycle) {
        Nextflow.error("ERROR: You must have EITHER an input_sample parameter OR an input_cycle parameter, but not both!")
    } else if(!params.input_sample && !params.input_cycle) {
        Nextflow.error("ERROR: You must have EITHER an input_sample parameter OR and input_cycle paramter!")
    }

    ch_versions = Channel.empty()

    ch_from_samplesheet.ashlar.view { "ashlar $it" }
    
    marker_sheet_index_map = make_marker_sheet_index_map(params.marker_sheet)
    ch_from_marker_sheet = Channel.fromSamplesheet("marker_sheet")
    //    .map { validate_marker_sheet(it, sample_sheet_index_map, marker_sheet_index_map) }

    // Format input for BASICPY
    // data_path = ch_from_samplesheet
    //     .map(it->"${it[1]}/*.ome.tif")
    // raw_cycles = Channel.of([[id:"exemplar-001"],"/workspace/data/exemplar-001/raw/exemplar-001-cycle-06.ome.tiff"])

    //
    // MODULE: BASICPY
    //

    // BASICPY(raw_cycles)
    //ch_versions = ch_versions.mix(BASICPY.out.versions)

    // /*
    // if ( params.illumination ) {
    //     BASICPY(ch_images)
    //     ch_tif = BASICPY.out.fields
    //     ch_versions = ch_versions.mix(BASICPY.out.versions)

    //     ch_dfp = ch_tif.filter { file -> file.name.endsWith('.dfp.tiff') }
    //     ch_ffp = ch_tif.filter { file -> file.name.endsWith('.ffp.tiff') }
    // }
    // */

    // MARKER_SHEET_CHECK(params.marker_sheet)

    // INPUT_CHECK(params.input_cycle, params.marker_sheet)
    INPUT_CHECK( input_type, params.input_sample, params.input_cycle, params.marker_sheet )
    // MARKER_CHECK(parmas.marker_sheet)

    // ASHLAR(raw_images, dfp, ffp)
    //ASHLAR(raw_images, [], [])
    ASHLAR(ch_from_samplesheet.ashlar, [], [])
    ch_versions = ch_versions.mix(ASHLAR.out.versions)

    // // Run Background Correction
    // BACKSUB(ASHLAR.out.tif, ch_markers)
    // ch_versions = ch_versions.mix(BACKSUB.out.versions)

    // // Run Segmentation
    /*
    DEEPCELL_MESMER(ASHLAR.out.tif, [[:],[]])
    ch_versions = ch_versions.mix(DEEPCELL_MESMER.out.versions)
    */

    /* starts running, but errors out. See ticket for error message
    def project = '/Users/robertyoung/Projects/test_mcmicro/Tonsil.ilp'
    ILASTIK_PIXELCLASSIFICATION( ASHLAR.out.tif, [[id:'test2'], project] )
    */

    // // Run Quantification
    /*
    MCQUANT(ASHLAR.out.tif,
            DEEPCELL_MESMER.out.mask,
            markerFile)
    ch_versions = ch_versions.mix(MCQUANT.out.versions)
    */
    // // Run Reporting
    // SCIMAP_MCMICRO(MCQUANT.out.csv)
    // ch_versions = ch_versions.mix(SCIMAP_MCMICRO.out.versions)

    /*
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    */

    //
    // MODULE: MultiQC
    //
    /*
    workflow_summary    = WorkflowMcmicro.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMcmicro.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
    */
}

def make_sample_sheet_index_map(String sample_sheet_path) {
    def sample_sheet_index_map = [:]
    def header
    new File(sample_sheet_path).withReader { header_list = it.readLine().split(',') }
    def ctr = 0
    header_list.each { value ->
        sample_sheet_index_map[value] = ctr
        ctr = ctr + 1
    }
    return sample_sheet_index_map
}

def make_marker_sheet_index_map(String marker_sheet_path) {
    def marker_sheet_index_map = [:]
    def header
    new File(marker_sheet_path).withReader { header_list = it.readLine().split(',') }
    def ctr = 0
    header_list.each { value ->
        marker_sheet_index_map[value] = ctr
        ctr = ctr + 1
    }
    return marker_sheet_index_map
}

def make_ashlar_input_sample(ArrayList sample_sheet_row, Map sample_sheet_index_map) {
    sample_name_index = sample_sheet_index_map['sample']
    image_dir_path_index = sample_sheet_index_map['image_directory']

    files = []
    def image_dir = new File(sample_sheet_row[image_dir_path_index])
    image_dir.eachFileRecurse (FileType.FILES) {
        if(it.toString().endsWith(".ome.tif")){
            files << file(it)
        }
    }

    ashlar_input = [[id:sample_sheet_row[sample_name_index]], files]

    return ashlar_input
}

def make_ashlar_input_cycle(ArrayList sample_sheet_row, Map sample_sheet_index_map) {
    sample_name_index = sample_sheet_index_map['sample']
    image_tiles_path_index = sample_sheet_index_map['image_tiles']
    ashlar_input = [[id:sample_sheet_row[sample_name_index]], sample_sheet_row[image_tiles_path_index]]

    return ashlar_input
}

marker_name_list = []
cycle_channel_tuple_list = []

/* moved validation to subworkflow
def validate_marker_sheet(ArrayList marker_sheet_row, Map sample_sheet_index_map, Map marker_sheet_index_map) {
    channel_index = marker_sheet_index_map['channel_number']
    cycle_index = marker_sheet_index_map['cycle_number']
    marker_index = marker_sheet_index_map['marker_name']
    sample_name_index = sample_sheet_index_map['sample']
    cycle_number_index = sample_sheet_index_map['cycle_number']
    channel_count_index = sample_sheet_index_map['channel_count']

    if (!sample_sheet_index_map['cycle_number']) {
        // we're currently only validating 1 row per sample per cycle sample sheets
        return
    }

    // check marker name uniqueness
    if(marker_name_list.contains(marker_sheet_row[marker_index])){
        Nextflow.error("ERROR: Duplicate marker name in marker sheet! Marker names must be unique.")
    } else {
        marker_name_list.add(marker_sheet_row[marker_index])
    }

    // check that (cycle, channel) tuple is unique
    curr_tuple = new Tuple(marker_sheet_row[channel_index], marker_sheet_row[cycle_index])
    if(cycle_channel_tuple_list.contains(curr_tuple)){
        Nextflow.error("ERROR: Duplicate cycle_number & channel_number pair! cycle_number & channel_number pairs must be unique.")
    } else if(marker_sheet_row[channel_index] < 1 || marker_sheet_row[cycle_index] < 1) {
        Nextflow.error("ERROR: cycle_number and channel_number are 1-based.  Values less than 1 are not allowed.")
    } else if((cycle_channel_tuple_list.size() > 0) &&
              ( ((marker_sheet_row[channel_index] != cycle_channel_tuple_list.last()[channel_index]) &&
                 (marker_sheet_row[channel_index] != cycle_channel_tuple_list.last()[channel_index]+1)) ||
                ((marker_sheet_row[cycle_index] != cycle_channel_tuple_list.last()[cycle_index]) &&
                 (marker_sheet_row[cycle_index] != cycle_channel_tuple_list.last()[cycle_index]+1)) )) {
        Nextflow.error("ERROR: cycle_number and channel_number must be sequential with no gaps")
    } else {
        cycle_channel_tuple_list.add(curr_tuple)
    }
}
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
