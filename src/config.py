
# General keys
EXT_DCM     = '.dcm'
MODALITY_CT = 'CT'

PREFIX_PATIENT_ID  = 'HCAI-Dose-{}'
PREFIX2_PATIENT_ID = 'HCAI-Dose-x{}'
PREFIX3_PATIENT_ID = 'HCAI-Dose-Pr{}'
PREFIX_PROTON_PATIENT_ID = 'HCAI-Dose-P{}'
PREFIX_PROTON2_PATIENT_ID = 'HCAI-Dose-Pr{}'

KEYNAME_CT       = 'CT'
KEYNAME_RTPLAN   = 'RTPLAN'
KEYNAME_RTDOSE   = 'RTDOSE'
KEYNAME_RTSTRUCT = 'RTSTRUCT'
KEYNAME_UNAPPROVED = 'UNAPPROVED'

KEYNAME_PATH_PATIENT    = 'pathPatient'
KEYNAME_PATH_CLASSSOL   = 'pathClassSolution'
KEYNAME_PATH_OBJECTIVES = 'pathKNOObjectivesClinical'
KEYNAME_PATH_DVHPARAMS  = 'pathDVHParams'
KEYNAME_PATH_ISODOSEXML = 'pathIsoDoseXML'
KEYNAME_OPT_STEPS_RE    = 'optStepsForRe'

KEYNAME_FORCE_LOAD_PATIENT    = 'forceLoadPatient'
KEYNAME_FORCE_UPLOAD_PATIENT  = 'forceUploadPatient'
KEYNAME_FORCE_CURRENT_PATIENT = 'forceCurrentPatient'

# RS-specific keys
KEYNAME_RS_PATIENTDB = 'PatientDB'
KEYNAME_PATIENT      = 'Patient'
KEYNAME_CASE         = 'Case'
KEYNAME_PTV          = 'PTV'
KEYNAME_RS_PLAN      = 'Plan'
KEYNAME_RS_BEAMSET   = 'BeamSet'
KEY_PATIENTID        = 'PatientID'

# LUMC-Specific keys
NAME_CT_IMAGING_SYSTEM = 'BBBCT120kV' # specific to LUMC
PHYSICIAN_NAME         = 'Mody' # comes up in the list of patients so that others researchers can identify
NAME_TREATMENT_MACHINE = 'Toestel 5' # specific to LUMC

SUFFIX_PLAN_CS     = '-{}1'
SUFFIX_PLAN_DFO    = '-{}2'
SUFFIX_PLAN_DFO2   = '-{}3'
SUFFIX_PLAN_EUD    = '-{}4'
SUFFIX_PLAN_FINAL  = '-{}5'
SUFFIX_PLAN_FINAL2 = '-{}6'

PREFIX_CLINICAL_CONTOURS  = 'R'
PREFIX_AUTOMATED_CONTOURS = 'A'

KEYNAME_CANCER_TYPE = 'cancerType'

KEYNAME_PLAN_OG     = '{}'
KEYNAME_PLAN_CS     = '{}' + SUFFIX_PLAN_CS # CLASS SOLUTION
KEYNAME_PLAN_DFO    = '{}' + SUFFIX_PLAN_DFO
KEYNAME_PLAN_DFO2   = '{}' + SUFFIX_PLAN_DFO2
KEYNAME_PLAN_EUD    = '{}' + SUFFIX_PLAN_EUD
KEYNAME_PLAN_FINAL  = '{}' + SUFFIX_PLAN_FINAL
KEYNAME_PLAN_FINAL2 = '{}' + SUFFIX_PLAN_FINAL2

KEYNAME_CLASS_SOLUTION_CONTOURS_OG   = 'ClsSol-OGCon'
KEYNAME_CLASS_SOLUTION_CONTOURS_AUTO = 'ClsSol-AutoCon'

KEYNAME_TREATMENT_MACHINE = 'Toestel 1'
KEYNAME_RADIATION_PHOTONS = 'Photons'
KEYNAME_TREATMENT_TECHNIQUE = 'VMAT'
KEYNAME_PLAN_COMMENT = "DL1 5245\nDL2 7000"

KEYNAME_BEAM_1 = "01.01"
KEYNAME_BEAM_1_DESC = "1 arc 178-182"
KEYNAME_BEAM_2 = "01.02"
KEYNAME_BEAM_2_DESC = "1 arc 182-178"

KEYNAME_ISO_CENTER = 'IsoCenter'

RS_CHECK_PREFIX = '_'
RS_STR_TRANSLATE_OBJ = {
    ord(' ') : ord('_')
    , ord('(') : ord('_')
    , ord(')') : ord('_')
    , ord('-') : ord('_')
    , ord('+') : ord('_')
    , ord('<') : ord('_')
    , ord('>') : ord('_')
}

KEY_AUTOCONTOUR_SUFFIX = ' (1)' # very rs-specific
KEY_3MM_SUFFIX = '+3'
KEY_OBJ_SUFFIX = '_obj'


KEYNAME_PTV_DL1_DVH = 'PTV_DL1_DVH'
KEYNAME_PTV_DL2_DVH = 'PTV_DL2_DVH'
KEYNAME_CAVITY_ORAL = 'Cavity_Oral' # RS Autocontouring uses AAPM nomenclature
KEYNAME_ORAL_CAVITY = 'Oral_Cavity' # Inside of KNO.xml
KEYNAME_BRAINSTEM = 'Brainstem'
KEYNAME_SPINALCORD = 'SpinalCord'
KEYNAME_PAROTID_L = 'Parotid_L'
KEYNAME_PAROTID_R = 'Parotid_R'
KEYNAME_GLND_SUBMAND_L = 'Glnd_Submand_L'
KEYNAME_GLND_SUBMAND_R = 'Glnd_Submand_R'
KEYNAME_SUBMAND_L_OBJ = 'Submand_L_obj'
KEYNAME_SUBMAND_R_OBJ = 'Submand_R_obj'
KEYNAME_CRICO = 'Cricopharyngeus'
KEYNAME_LARYNX_SG = 'Larynx_SG'
KEYNAME_GLOTTIC_AREA = 'Glottic_Area'
KEYNAME_COCHLEA_L = 'Cochlea_L'
KEYNAME_COCHLEA_R = 'Cochlea_R'

KEYNAME_MUSC_CONSTRICT_I = 'Musc_Constrict_I'
KEYNAME_MUSC_CONSTRICT_M = 'Musc_Constrict_M'
KEYNAME_MUSC_CONSTRICT_S = 'Musc_Constrict_S'
KEYNAME_ESOPHAGUS        = 'Esophagus' # <-- NEED TO RESOLVE THIS VISUALLY
KEYNAME_ESOPHAGUS_S      = 'Esophagus_S' # <-- NEED TO RESOLVE THIS VISUALLY
KEYNAME_TRACHEA          = 'Trachea'
KEYNAME_SWAL_COMP        = 'Swal_Comp' # KEYNAME_MUSC_CONSTRICT_I + KEYNAME_MUSC_CONSTRICT_M + KEYNAME_MUSC_CONSTRICT_S + KEYNAME_CRICO + KEYNAME_LARYNX_SG + KEYNAME_GLOTTIC_AREA
KEYNAME_SWAL_OBJ         = 'Swal_obj'

KEYNAME_MANDIBLE = 'Bone_Mandible'
KEYNAME_MANDIBLE_PTV = 'Bone_Mandible-PTV'
KEYNAME_RING_LT_PTV_DL1 = 'ring<PTV_DL1'
KEYNAME_BRAIN = 'Brain'


KEYNAME_GHOST_CRANIAL = 'Ghost_cranial'
KEYNAME_EAR_L_GHOST = 'Ear_L_ghost'
KEYNAME_EAR_R_GHOST = 'Ear_R_ghost'

KEYNAME_BODY         = 'Body'
KEYNAME_OPT_BODY     = 'Opt_Body'
KEYNAME_RING_PTV_DL2 = 'ring>PTV_DL2'
KEYNAME_GHOST        = 'ghost'

REGEX_PROSTHESE    = 'prosthese'
REGEX_DFO          = 'dfo'
REGEX_DMAX         = 'dmax'
REGEX_DOSE         = 'dose'
REGEX_KAAK         = 'kaak'
REGEX_LIPPEN       = 'lippen'
REGEX_5805         = '5805'



KEY_XML_ROOT_NAME = 'RayStationScriptTemplates'
KEY_XML_SUBROOT_NAME = 'ObjectiveTemplate'
KEY_XML_ROI = 'Roi'
KEY_XML_NAME = 'name'

KEY_FTYPE_DOSEFALLOFF = 'DoseFallOff'
KEY_FTYPE_MINDOSE = 'MinDose'
KEY_FTYPE_MAXDOSE = 'MaxDose'
KEY_FTYPE_UNIFORMDOSE = 'UniformDose'
KEY_FTYPE_MINDVH = 'MinDvh'
KEY_FTYPE_MAXDVH = 'MaxDvh'
KEY_FTYPE_MAXEUD = 'MaxEud'
KEY_FTYPE_MINEUD = 'MinEud'
KEY_FTYPE_UNIFORMEUD = 'UniformEud'
KEY_FTYPE_UNIFORMITYCONSTRAINT = 'UniformityConstraint'

# OARS =  [KEYNAME_BRAINSTEM, KEYNAME_COCHLEA_L, KEYNAME_COCHLEA_R, KEYNAME_ESOPHAGUS_S, "Eye_L", "Eye_R", "Glnd_Lacrimal_L", "Glnd_Lacrimal_R", "Glottis"
#                     , KEYNAME_LARYNX_SG, "Lens_L", "Lens_R",  "Bone_Mandible",  "Nasolacrimal_Duct_L", "Nasolacrimal_Duct_R", "Nasopharynx", "OpticNrv_L", "OpticNrv_R", KEYNAME_CAVITY_ORAL, "Oropharynx"
#                     , KEYNAME_PAROTID_L, KEYNAME_PAROTID_R, "Pituitary", "Fossa_Posterior", KEYNAME_SPINALCORD, KEYNAME_GLND_SUBMAND_L, KEYNAME_GLND_SUBMAND_R, "Joint_TM_L", "Joint_TM_R", "Tongue_Base"]

# This is a union of items from the KNO.xml and the RS Autocontouring Module
OARS =  [KEYNAME_BRAINSTEM, KEYNAME_SPINALCORD
         , KEYNAME_COCHLEA_L, KEYNAME_COCHLEA_R
         , KEYNAME_ESOPHAGUS_S, KEYNAME_LARYNX_SG, KEYNAME_CAVITY_ORAL, KEYNAME_MANDIBLE
         , KEYNAME_PAROTID_L, KEYNAME_PAROTID_R, KEYNAME_GLND_SUBMAND_L, KEYNAME_GLND_SUBMAND_R
        ] 
OAR_DUPLICATE_COLOR_RGB_STRING = "128,0,64" # purple
OAROBJ_DUPLICATE_COLOR_RGB_STRING = "255,128,128" # pink
# [Frank: 2023-06-02] in the clinic: Brainstem, Esophagus,Bone_Mandible,Oral_Cavity,SpinalCord, Glnd_{Par/SMD}_{L/R}

PHOTON_POTENTIAL_ROIS_TO_RENAME_FOR_AUTO = [
        KEYNAME_BRAINSTEM, KEYNAME_BRAINSTEM + KEY_3MM_SUFFIX
        , KEYNAME_SPINALCORD, KEYNAME_SPINALCORD + KEY_3MM_SUFFIX
        , KEYNAME_PAROTID_L, KEYNAME_PAROTID_L + KEY_OBJ_SUFFIX
        , KEYNAME_PAROTID_R, KEYNAME_PAROTID_R + KEY_OBJ_SUFFIX
        , KEYNAME_ORAL_CAVITY, KEYNAME_ORAL_CAVITY + KEY_OBJ_SUFFIX
        , KEYNAME_GLND_SUBMAND_L, KEYNAME_SUBMAND_L_OBJ
        , KEYNAME_GLND_SUBMAND_R, KEYNAME_SUBMAND_R_OBJ
        , KEYNAME_COCHLEA_L, KEYNAME_COCHLEA_R
        , KEYNAME_LARYNX_SG
        , KEYNAME_SWAL_OBJ
        , KEYNAME_SWAL_COMP
        , KEYNAME_ESOPHAGUS
        , KEYNAME_MANDIBLE, KEYNAME_MANDIBLE_PTV
    ]

###########################################################################
# PROTON PLANNING
###########################################################################

KEYNAME_ROI_CTV_DL1               = 'CTV_DL1'
KEYNAME_3MM                       = '3mm'
KEYNAME_ROI_CTV_DL1_3MM           = '{}+{}'.format(KEYNAME_ROI_CTV_DL1, KEYNAME_3MM)
KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX = '-({})'.format(KEYNAME_ROI_CTV_DL1_3MM) #
KEYNAME_ROI_MID_STRUCTURES        = 'Mid_structures'

REGEX_ROI_HOT       = 'hot'
REGEX_ROI_HEET      = 'heet'
REGEX_ROI_HOSTPSOTS = 'hotspot' 
REGEX_ROI_MINDOSE   = 'min dose'
REGEX_ROI_MAXDOSE   = 'max dose'
REGEX_DFO           = 'dfo'
REGEX_ONDER         = 'onder'
REGEX_OVER          = 'over'
REGEX_KOUD          = 'koud'
REGEX_COLD          = 'cold'
REGEX_ONDERDOSERING = 'onderdosering'
REGEX_CAUDAAL       = 'caudaal'
REGEX_CAUD          = 'caud'

KEYNAME_ROI_PTV_DL1 = 'PTV_DL1'
KEYNAME_ROI_PTV_DL2 = 'PTV_DL2'

PROTON_POTENTIAL_ROIS_TO_RENAME_FOR_AUTO = [
        KEYNAME_BRAINSTEM, KEYNAME_BRAINSTEM + KEY_3MM_SUFFIX
        , KEYNAME_SPINALCORD, KEYNAME_SPINALCORD + KEY_3MM_SUFFIX
        , KEYNAME_PAROTID_L, KEYNAME_PAROTID_L + KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX
        , KEYNAME_PAROTID_R, KEYNAME_PAROTID_R + KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX
        , KEYNAME_GLND_SUBMAND_L, KEYNAME_GLND_SUBMAND_L + KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX
        , KEYNAME_GLND_SUBMAND_R, KEYNAME_GLND_SUBMAND_R + KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX
        , KEYNAME_ORAL_CAVITY, KEYNAME_ORAL_CAVITY + KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX
        , KEYNAME_COCHLEA_L, KEYNAME_COCHLEA_R
        , KEYNAME_LARYNX_SG
        , KEYNAME_ROI_MID_STRUCTURES + KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX
        , KEYNAME_ESOPHAGUS
        , KEYNAME_MANDIBLE
    ]

FILENAME_ROBUST_EVAL_RESULTS = 'robustEvalResults-{}.json'
KEYNAME_PATH_ROBUST_TEMPLATE = 'pathRobustEvalTemplate'

KEYNAME_PASS    = 'Pass'
KEYNAME_FAIL    = 'Fail'
KEYNAME_ATMOST  = 'AtMost'
KEYNAME_ATLEAST = 'AtLeast'

CTV_DL1_MAX = 5425
CTV_DL2_MAX = 7000

KEYNAME_SCENARIO_DOSE = 'ScenarioDose'
KEYNAME_NOMINAL_DOSE  = 'NominalDose'
KEYNAME_MAX_DOSE      = 'MaxDose'
KEYNAME_MIN_DOSE      = 'MinDose'
KEYNAME_PASSED        = 'Passed'
KEYNAME_VOXELWISE_WORST = 'Voxelwise-worst'
KEYNAME_SCENARIOS       = 'Scenarios'
OBJECTIVES_ROIS_TUMOR = ['CTV_DL1 [D98%>95%]', 'CTV_DL1 [D98%>94%]', 'CTV_DL2 [D98%>95%]', 'CTV_DL2 [D98%>94%]']


###########################################################################
# PLAN STATS
###########################################################################
KEYNAME_PLAN_DEBUGINFO   = 'planDebugInfo'
KEYNAME_TIME             = 'time'
KEYNAME_OBJ_VALUE        = 'objValue'
KEYNAME_PLAN_EXTRAS      = 'planExtras'
KEYNAME_AUTOCONTOURING_TIME = 'autoContouringTime'
KEYNAME_PLANS            = 'plans'
KEY_PLAN_STATS_JSON      = 'planStats.json'
KEY_PLAN_STATS_JSON_AUTO = 'planStatsAutoContour.json'
KEY_PLAN_STATS_JSON_ALL  = 'planStatsAll.json'
KEYNAME_CONTOUR_TYPE     = 'contourType'
KEYNAME_CONTOUR_CLINICAL = 'clinical-contour'
KEYNAME_CONTOUR_AUTO     = 'auto-contour'
KEYNAME_CONTOUR_ALL      = 'all-contour'
KEYNAME_CONTOUR_EVAL     = 'eval-contour'
KEYNAME_CONTOUR_DEBUG     = 'eval-debug'

######################################################################
# NTCP KEYS 
######################################################################

KEY_NTCP_PLAN1 = "plan_1"
KEY_NTCP_PLAN2 = "plan_2"
KEY_NTCP_TREATMENT_TYPE = "treatment_type"
KEY_NTCP_TUMOR_LOCATION = "tumor_location"
KEY_NTCP_BASELINE_XEROSTOMIA   = "baseline_xerostomia"
KEY_NTCP_BASELINE_DYSPHAGIA    = "baseline_dysphagia"
KEY_NTCP_IS_TOTAL_LARYNGECTOMY = "is_total_laryngectomy"
KEY_PAROTIDS_REMOVED           = "parotids_removed"

KEY_NTCP_PLAN_1 = "Plan_1"
KEY_NTCP_PLAN_VALUE = "PlanValue"
KEY_NTCP_MODEL_NAME = "ModelName"

# FILENAME_NTCP_RESULTS = 'ntcpResults.json'
FILENAME_NTCP_RESULTS = 'ntcpResultsV2.json'
FILENAME_NTCP_PHOTON_RESULTS = 'ntcpPhotonResults.json'

KEY_XERO_GRADE2 = 'Xerostomia Grade ≥ 2'
KEY_XERO_GRADE3 = 'Xerostomia Grade ≥ 3'
KEY_DYS_GRADE2  = 'Dysphagia Grade ≥ 2'
KEY_DYS_GRADE3  = 'Dysphagia Grade ≥ 3'


"""
CTV_DL2_DVH (required)
CTV_DL1_DVH (required)
PTV_DL2_DVH (required)
PTV_DL1_DVH (required)
Brainstem (required)
Brain
SpinalCord (required)
Parotid_L (required)
Parotid_R (required)
Glnd_Submand_L
Glnd_Submand_R
Musc_Constrict_I (required)
Musc_Constrict_M (required)
Musc_Constrict_S (required)
Cricopharyngeus (required)
Oral_Cavity
Esophagus (required)
Bone_Mandible
Bone_Mandible-PTV
Trachea (required)
Glottic_Area (required)
Larynx_SG (required)
PTV_DL1_obj_5mm (required)
PTV_DL1_obj_10mm (required)
Cochlea_L (required)
Cochlea_R (required)
Lens_L
Lens_R
Ghost_craniaal (required)
Body (required)
"""

"""
allowed additionals:
- planNameClassSol      : Ghost_craniaal, Ear_L_ghost, Ear_R_ghost
- planNameCSAndDFO      : N/A
- planNameCSAndDFOAndEUD: N/A
- planNameFinalTouches  : prothese, dfo, dmax, dose, ring>PTV_DL2
"""