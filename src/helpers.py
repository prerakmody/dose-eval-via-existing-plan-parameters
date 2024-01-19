# Import RStation libraries
import connect

# Import private modules
import hnDoseConfig as config

# Import public modules
import re
import os
import pdb
import copy
import time
import math
import shutil
import logging
import pydicom
import traceback
import numpy as np
from pathlib import Path
from xml.etree import ElementTree
from distutils.util import strtobool

def print(*args, **kwargs):
    logging.info(" ".join(map(str, args)), **kwargs)

##########################################################################################
#                                        MY CODE                                         #
##########################################################################################

########## PATIENT-RELATED ########## 

def getPatientIdentifier(patientObj):
    """
    patientObj: connect.connect_cpython.PyScriptObject
    """
    return '/'.join([patientObj.PatientID, patientObj.Name])

def rayStationSave():

    patientObj = None

    try:
        patientObj = connect.get_current(config.KEYNAME_PATIENT)
        patientObj.Save() # need to do this so to avoid "PreConditionViolationException: State must be saved."
    except:
        pass # if there is no patient open
    
    return patientObj

def getPatientById(patientID, lastFind=True):

    patient = None

    db        = connect.get_current(config.KEYNAME_RS_PATIENTDB)
    patients  = db.QueryPatientInfo(Filter={config.KEY_PATIENTID: patientID})
    patients  = sorted(patients, key=lambda x: x['LastName'], reverse=False)
    lastNames = [each['LastName'] for each in patients]
    if len(patients):
        if lastFind:
            print (f" - [getPatientById] Loading {patients[-1]['LastName']} from {lastNames}")
            patient = db.LoadPatient(PatientInfo=patients[-1], AllowPatientUpgrade=True)
        else:
            print (f" - [getPatientById] Loading {patients[0]['LastName']} from {lastNames} ")
            patient = db.LoadPatient(PatientInfo=patients[0], AllowPatientUpgrade=True)
    
    return patient

def loadPatientUsingID(patientId):
    """
    Params
    ------
    patientId: str
    """

    patientRSObj = None

    try:
        db = connect.get_current(config.KEYNAME_RS_PATIENTDB)
        patientInfos = db.QueryPatientInfo(Filter={config.KEY_PATIENTID: '^' + str(patientId) + '$'}) # regex filtering

        if len(patientInfos):
            # [NOTE: Returns the most recent patient]
            patientRSObj = db.LoadPatient(PatientInfo=patientInfos[-1], AllowPatientUpgrade=False)
    except:
        traceback.print_exc()

    return patientRSObj

########## DICOM-RELATED (upload) ##########

def uploadRTAppsDataToRStation(pathPatient, planName, forceUpload=False, forceCurrentPatient=False):
    """
    Params
    ------
    pathPatient: Path, Path to patient folder containing CT, RTDose, RTPlan, RTStruct
     - e.g. Path('H:\\').joinpath('RayStationData', 'LUMC-Dose', 'HCAI-Dose-1', '2.25.52093085334020578701550802539023326955')
        - CT_{}
        - RTDOSE_{}
        - RTPLAN_{}
        - RTSTRUCT_{}
    planName: String, check if this plan exists in RayStation. If not, upload
    forceUpload: Bool, If True, will upload data even if patient already exists. Useful when debugging
    forceCurrentPatient: Bool, If True, will use patient currently open in RayStation
    
    NOTE: Above folder structure is required for this function to work
    """
    
    # Step 0 - Init
    assert forceUpload + forceCurrentPatient < 2, ' - [uploadRTAppsDataToRStation] forceUpload and forceCurrentPatient cannot be both True!'
    rayStationSave()
    db = connect.get_current(config.KEYNAME_RS_PATIENTDB)
    patientCTBool, patientRTStructBool, patientRTPlanBool, patientRTDoseBool = False, False, False, False
    
    # Step 1 - Create patient
    if pathPatient.exists():

        # Step 2 - Get patient paths
        pathPatientCTFolders       = [path for path in pathPatient.iterdir() if path.is_dir() and path.name.startswith(config.KEYNAME_CT)]
        pathPatientRTDoseFolders   = [path for path in pathPatient.iterdir() if path.is_dir() and path.name.startswith(config.KEYNAME_RTDOSE)]
        pathPatientRTPlanFolders   = [path for path in pathPatient.iterdir() if path.is_dir() and path.name.startswith(config.KEYNAME_RTPLAN)]
        pathPatientRTStructFolders = [path for path in pathPatient.iterdir() if path.is_dir() and path.name.startswith(config.KEYNAME_RTSTRUCT)]

        if len(pathPatientCTFolders) == 1 and len(pathPatientRTDoseFolders) == 1 and len(pathPatientRTPlanFolders) == 1 and len(pathPatientRTStructFolders) == 1:
            pathPatientCTFolder       = pathPatientCTFolders[0]
            pathPatientRTDoseFolder   = pathPatientRTDoseFolders[0]
            pathPatientRTPlanFolder   = pathPatientRTPlanFolders[0]
            pathPatientRTStructFolder = pathPatientRTStructFolders[0]

            # Step 1.1 - Upload CT
            patientIdCheck = Path(pathPatient).parts[-2]
            if not forceCurrentPatient:
                try:
                    print ('\n\n - [uploadRTAppsDataToRStation()] --------------------- Step 1: Checking for patientId={}, plan={} (forceUpload={}) \n\n'.format(patientIdCheck, planName, forceUpload))
                    if forceUpload: patient = None
                    else          : patient = getPatientById(patientIdCheck, lastFind=True)
                    if patient is None:
                        print ('\n - [uploadRTAppsDataToRStation()] --------------------- Step 1.1: Uploading patient data for patientId={} \n\n'.format(patientIdCheck))
                        if Path(pathPatientCTFolder).exists():
                            patientID, studyUID, seriesUID = updateCTDicoms(pathPatientCTFolder)
                            print ('\n - [uploadRTAppsDataToRStation()] --------------------- Step 1.1: Uploading CT data for {} ... \n\n'.format(patientID))
                            warnings = db.ImportPatientFromPath(Path=str(pathPatientCTFolder), SeriesOrInstances=[{'PatientID': patientID, 'StudyInstanceUID': str(studyUID), 'SeriesInstanceUID': str(seriesUID)}], ImportFilter='', BrachyPlanImportOverrides={})
                            patient = rayStationSave()
                            # patient.Cases[0].SetCurrent()
                            print ('\n - [uploadRTAppsDataToRStation()] --------------------- Step 1.1: Uploaded CT data for {} \n\n'.format(getPatientIdentifier(patient)))
                            setEquipmentName()
                            patientCTBool = True
                        else:
                            print ('\n - [uploadRTAppsDataToRStation()][{}] --------------------- Step 1.1: No CT data found: {} \n\n'.format(patientID, pathPatientCTFolder))
                    else:
                        print ('\n - [uploadRTAppsDataToRStation()] --------------------- Step 1.1: Patient already exists \n\n')
                        patientCTBool = True
                except:
                    traceback.print_exc()
                    print ('\n - [uploadRTAppsDataToRStation()] --------------------- Step 1.1: Issue with CT data upload: {} \n\n'.format(pathPatientCTFolder))
                    patientCTBool = False
            else:
                print ('\n\n [uploadRTAppsDataToRStation()][DEBUG] --------------------- Step 1: Not uploading for patientID={} (forceCurrentPatient={}) \n\n'.format(patientIdCheck, forceCurrentPatient))
                patientCTBool = True
            print ('\n ---------------------------------------------------------- ')

            # Step 1.2 - Upload (existing) RTStruct/RTPLAN/RTDOSE
            patient = rayStationSave()
            if patient is not None and patientCTBool:
                case      = patient.Cases[0]
                casename  = case.CaseName
                patientID = patient.PatientID
                
                # Step 1.2 - Upload (existing) RTStruct
                try:
                    print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: Upload RTStruct data \n\n'.format(getPatientIdentifier(patient)))
                    if not checkForRTStruct(case):
                        pathPatientRTStructFile = [each for each in pathPatientRTStructFolder.iterdir()][0]
                        if Path(pathPatientRTStructFile).exists():
                            studyUIDRTStruct, seriesUIDRTStruct = updateRTStructDicoms(pathPatientRTStructFile)
                            print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: Uploading RTStruct data ... \n\n'.format(getPatientIdentifier(patient)))
                            warningsRTStruct  = patient.ImportDataFromPath(Path=str(pathPatientRTStructFolder), CaseName=casename, SeriesOrInstances=[{'PatientID': patientID, 'StudyInstanceUID': str(studyUIDRTStruct), 'SeriesInstanceUID': str(seriesUIDRTStruct)}], ImportFilter='', BrachyPlanImportOverrides={}, AllowMismatchingPatientID=True)
                            print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: Uploaded RTStruct data \n\n'.format(getPatientIdentifier(patient)))
                            patientRTStructBool = True
                            rayStationSave()
                        else:
                            print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: No RTStruct data found: {} \n\n'.format(getPatientIdentifier(patient), pathPatientRTStructFile))
                    else:
                        print ('\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: RTStruct data already exists \n\n'.format(getPatientIdentifier(patient)))
                        patientRTStructBool = True
                except:
                    traceback.print_exc()
                    print ('\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: Issue with RTStruct data upload: {} \n\n'.format(getPatientIdentifier(patient), pathPatientRTStructFile))
                    patientRTStructBool = False

                print ('\n ---------------------------------------------------------- ')
                    
                # Step 1.3 - Upload (existing) RTPlan/RTDose
                if patientCTBool and patientRTStructBool:
                    if not checkForRTPlan(case, planName):
                        print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.3: Upload RTPlan/RTDose data \n\n'.format(getPatientIdentifier(patient)))
                        studyUIDRTPlan  = Path(pathPatientRTPlanFolder).parts[-2]
                        seriesUIDRTPlan = Path(pathPatientRTPlanFolder).parts[-1].split('_')[-1]
                        studyUIDRTDose  = Path(pathPatientRTDoseFolder).parts[-2]
                        seriesUIDRTDose = Path(pathPatientRTDoseFolder).parts[-1].split('_')[-1]
                        updateRTPlanDicoms(pathPatientRTPlanFolder)
                        pathTempRTDoseAndRTPlanFolder = getTempRTDoseAndRTPlanFolder(pathPatientRTPlanFolder, pathPatientRTDoseFolder)
                        print ('\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.3: Uploading RTPlan/RTDose data ... \n\n'.format(getPatientIdentifier(patient)))
                        warningsRTPlan  = patient.ImportDataFromPath(Path=str(pathTempRTDoseAndRTPlanFolder), CaseName=casename, SeriesOrInstances=[
                            {'PatientID': patientID, 'StudyInstanceUID': str(studyUIDRTPlan), 'SeriesInstanceUID': str(seriesUIDRTPlan)}    
                            , {'PatientID': patientID, 'StudyInstanceUID': str(studyUIDRTDose), 'SeriesInstanceUID': str(seriesUIDRTDose)}
                            ], ImportFilter='', BrachyPlanImportOverrides={}, AllowMismatchingPatientID=True)
                        shutil.rmtree(pathTempRTDoseAndRTPlanFolder)
                        patientRTPlanBool = True
                        patientRTDoseBool = True
                        print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.3: Uploaded RTPlan/RTDose data \n\n'.format(getPatientIdentifier(patient)))
                    else:
                        print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.3: RTPlan/RTDose already exists \n\n'.format(getPatientIdentifier(patient)))
                        patientRTPlanBool = True
                        patientRTDoseBool = True
                else:
                    print ('\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.3: Not uploading RTPlan/RTDose due to CT/RTStruct issue \n\n'.format(patientID))
                print ('\n ---------------------------------------------------------- ')

            else:
                if patientCTBool is False:
                    print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: No CT to be found for: \n\n'.format(pathPatientCTFolder))
                else:
                    print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: No patient to be found: \n\n'.format(patientID))
        
        else:
            print(' - [uploadRTAppsDataToRStation()] Patient folder does not contain one CT, RTDOSE, RTPLAN and RTSTRUCT folder: ', pathPatient)
    
    else:
        print(" - [uploadRTAppsDataToRStation()] Patient folder does not exist: ", pathPatient)
    
    if patientCTBool and patientRTStructBool and patientRTPlanBool and patientRTDoseBool:
        return True
    else:
        return None

########## DICOM-RELATED ##########

def updateRTStructDicoms(pathRTStruct):
    
    # Step 1 - Get the study and series UID
    dsRTStruct        = pydicom.dcmread(str(pathRTStruct))
    studyUIDRTStruct = dsRTStruct.StudyInstanceUID
    seriesUIDRTStruct = dsRTStruct.SeriesInstanceUID
    dsRTStruct.ApprovalStatus = config.KEYNAME_UNAPPROVED # [UNAPPROVED, APPROVED, REJECTED]
    
    # Step 2 - Update the RTROIObservationsSequence
    for structureSet in dsRTStruct.StructureSetROISequence:
        if structureSet.ROIName in [config.KEYNAME_PTV_DL1_DVH, config.KEYNAME_PTV_DL2_DVH]:
            roiNumber = structureSet.ROINumber
            print (roiNumber, structureSet.ROIName)
            for roiObservation in dsRTStruct.RTROIObservationsSequence:
                if roiObservation.ReferencedROINumber == roiNumber:
                    roiObservation.RTROIInterpretedType = config.KEYNAME_PTV # [PTV, CTV, GTV, ORGAN, ]
            # dsRTStruct.RTROIObservationsSequence[roiNumber-1].RTROIInterpretedType = config.KEYNAME_PTV # [PTV, CTV, GTV, ORGAN, ]

    # Step 3 - Save the file
    dsRTStruct.save_as(str(pathRTStruct), write_like_original=True)
    
    return studyUIDRTStruct, seriesUIDRTStruct

def updateCTDicoms(pathCT):

    patientID, studyUID, seriesUID = None, None, None

    for pathDicomCT in pathCT.iterdir():
        if pathDicomCT.suffix == config.EXT_DCM:
            
            ds = pydicom.dcmread(str(pathDicomCT), force=True)

            if ds.Modality == 'CT':
                patientID = ds.PatientID
                studyUID = ds.StudyInstanceUID
                seriesUID = ds.SeriesInstanceUID

                ds.ReferringPhysicianName = config.PHYSICIAN_NAME # NOTE: This is the update

                ds.save_as(str(pathDicomCT), write_like_original=True)
    
    return patientID, studyUID, seriesUID

def updateRTPlanDicoms(pathPatientRTPlanFolder):

    studyID = None

    if Path(pathPatientRTPlanFolder).exists():
        pathPatientRTPlanFile = [each for each in pathPatientRTPlanFolder.iterdir()][0]
        ds = pydicom.filereader.dcmread(str(pathPatientRTPlanFile))

        ds.ApprovalStatus = config.KEYNAME_UNAPPROVED # [UNAPPROVED, APPROVED, REJECTED]
        print (' - [updateRTPlanDicoms()] Updating beams to : ', config.NAME_TREATMENT_MACHINE)
        for beam in ds.BeamSequence:
            beam.TreatmentMachineName = config.NAME_TREATMENT_MACHINE
        
        studyID = str(ds.StudyInstanceUID)

        ds.save_as(str(pathPatientRTPlanFile), write_like_original=True)

    else:
        raise Exception(f" - [updateRTPlanDicoms] {pathPatientRTPlanFolder} does not exist!")

    return studyID

def setEquipmentName():
    examination = connect.get_current('Examination')
    examination.EquipmentInfo.SetImagingSystemReference(ImagingSystemName=config.NAME_CT_IMAGING_SYSTEM)
    rayStationSave()

def checkForRTStruct(case):
    exists = '_' + config.KEYNAME_PTV_DL1_DVH in dir(case.PatientModel.StructureSets[0].RoiGeometries)
    return exists

########## RS PLAN RELATED ##########

def getPatientAndPlan(planName, debug=False):
    
    patient, case, plan, beamset = None, None, None, None
    
    # Step 1 - Get patient
    patient = rayStationSave()
    if patient is None:
        return patient, case, plan, beamset
    case = patient.Cases[0]
    
    # Step 2 - Get plan
    if not checkForRTPlan(case, planName):
        return patient, case, plan, beamset
    plan    = case.TreatmentPlans[planName]
    beamset = plan.BeamSets[planName]

    # Step 3 - Return
    return patient, case, plan, beamset

def copyPlan(basePlanName, newPlanName, createArcBeam=True, debug=False):

    # Step 1 - Init
    copyPlanStatus = False

    # Step 2 - Boilerplate code for get patient and plans
    patient, case, basePlan, _ = getPatientAndPlan(basePlanName)
    patient, case, newPlan, _  = getPatientAndPlan(newPlanName)

    if patient is None:
        print (' - [ERROR][copyPlan()] No patient loaded')
        return copyPlanStatus
    
    if basePlan is None:
        print (' - [ERROR][copyPlan()] No basePlan named {0} found'.format(basePlanName))
        return copyPlanStatus
    
    # Step 3 - Copy plan
    if newPlan is None:
        
        print (' - [INFO][copyPlan()] Copying plan {0} to {1}'.format(basePlanName, newPlanName))
        case.CopyPlan(PlanName=basePlanName, NewPlanName=newPlanName, KeepBeamSetNames=False)
        rayStationSave()
        print (' - [INFO][copyPlan()] Copied plan: ', newPlanName)

        if createArcBeam:
            arcBeamStatus = makeDualArcBeam(newPlanName)
            if not arcBeamStatus:
                print (' - [ERROR][copyPlan()] Failed to make dual arc beam')
                return copyPlanStatus
            rayStationSave()
        
    else:
        print (' - [WARNING][copyPlan()] Plan named {0} already exists'.format(newPlanName))
    
    # Step 4 - Debug
    if debug:
        print (' - [INFO][copyPlan()] Setting current plan: ', newPlanName)
        patient.Cases[0].SetCurrent()
        patient.Cases[0].TreatmentPlans[newPlanName].SetCurrent()

    copyPlanStatus = True
    return copyPlanStatus

def getTempRTDoseAndRTPlanFolder(pathPatientRTPlan, pathPatientRTDose):
    """
    We have to copy the RTDose and RTPlan files to a temporary folder because the RayStation API does not allow an APPROVED RTPlan to have a RTDose uploaded to it!
    """
    
    import shutil
    import tempfile
    
    # Step 1 - Get files to copy
    pathPatientRTPlanFile = [each for each in pathPatientRTPlan.iterdir()][0]
    pathPatientRTDoseFile = [each for each in pathPatientRTDose.iterdir()][0]

    tmpDir = Path(pathPatientRTDose).joinpath('tmp')
    if not Path(tmpDir).exists():
        Path(tmpDir).mkdir(parents=True, exist_ok=True)
        shutil.copy(str(pathPatientRTPlanFile), str(tmpDir))
        shutil.copy(str(pathPatientRTDoseFile), str(tmpDir))

    return tmpDir

def checkForRTPlan(case, planName):
    
    # Step 0 - Init
    exists = False

    # Step 1 - Convert to RS naming standard
    planNameRSStandard = config.RS_CHECK_PREFIX + str(planName).translate(config.RS_STR_TRANSLATE_OBJ)
    if planNameRSStandard in dir(case.TreatmentPlans):
        exists = True
    
    return exists

def getRTPlanIndex(case, planName):
    
    # Step 0 - Init
    index = -1

    # Step 1 - Convert to RS naming standard
    try:
        index     = [plan.Name for plan in case.TreatmentPlans].index(planName)
    except:
        pass

    return index

def checkForRTDose(beamset):
    exists = False
    if beamset.FractionDose.DoseValues is not None:
        if len(beamset.FractionDose.DoseValues.AlgorithmProperties.DoseAlgorithm) > 0:
            exists = True
    
    return exists

def applyIsoDoseColors(pathIsoDoseXML, case, beamset):

    try:
        if pathIsoDoseXML is not None:
            if Path(pathIsoDoseXML).exists():
                doselevels = {dl: int(value) for dl, value in re.findall("(DL[0-9]) ([0-9]+)", beamset.Comment)}
                isoman = IsodoseManager(pathIsoDoseXML, case, doselevels)
                isoman.map_isodose()
            else:
                print (' - [copyPlanAndReRunOptimization()] --------------------- Step 3: IsoDose XML does not exist: ', pathIsoDoseXML)
    except:
        traceback.print_exc()

def doPlanSanityCheckOld(case, planName):
    
    """
    1. Remove rois in objective list that do not have a voume
    """

    # Step 0 - Init
    planObj = case.TreatmentPlans[planName]

    # Step 1 - Loop over all rois
    for roi in case.PatientModel.RegionsOfInterest:
        roiName = roi.Name
        volume = getRoiVolume(case, roiName)

        # Step 2 - If volume is less than 1e-6, remove from objective list of plan
        if volume < 1e-6:
            
            # Remove from plan objectives
            for roiFunc in planObj.PlanOptimizations[0].Objective.ConstituentFunctions:
                if roiFunc.ForRegionOfInterest.Name == roiName:
                    print ('\n - [doPlanSanityCheck()][{}] Removing roi: {} (vol={:.3f}) from objective list: '.format(planName, roi.Name, volume))
                    roiFunc.DeleteFunction()
                    print ('   - Deleted!')

def doPlanSanityCheck(case, planName):
    
    """
    1. Remove rois in objective list that do not have a voume
    """

    # Step 0 - Init
    planObj = case.TreatmentPlans[planName]

    # Step 1 - Loop over all rois (in objectives)
    for roiFunc in planObj.PlanOptimizations[0].Objective.ConstituentFunctions:
        roiName = roiFunc.ForRegionOfInterest.Name
        volume = getRoiVolume(case, roiName)

        # Step 2 - If volume is less than 1e-6, remove from objective list of plan
        if volume < 1e-6:
                
            # Remove from plan objectives
            print ('  - [doPlanSanityCheck()][{}] Removing roi function from objective list: {} (vol={:.3f}) from objective list: '.format(planName, roiName, volume))
            roiFunc.DeleteFunction()
            print ('   - Deleted!')

########## BEAM RELATED ##########

def makeDualArcBeam(planName):

    print (' - [INFO][makeDualArcBeam()] Making dual arc beam for plan: ', planName)
    try:
        # Step 1 - Get basics
        patient = connect.get_current(config.KEYNAME_PATIENT)
        plan = patient.Cases[0].TreatmentPlans[planName]
        beamset = plan.BeamSets[planName]

        # Step 2 - Delete all existing beams (NOTE: Does that also delete dose ??)
        isoCenter = None
        while len(beamset.Beams) != 0:
            for beam in beamset.Beams:
                isoCenter = beam.Isocenter.Position                
                beamset.DeleteBeam(BeamName=beam.Name)
        print (' - [INFO][makeDualArcBeam()] isoCenter: ', isoCenter)

        # Step 3 - Create new beam
        if isoCenter is None:
            return False
        
        beamset.CreateDefaultIsocenterData(Position=isoCenter)
        beam = beamset.CreateArcBeam(Name="01.01", Description="1 arc 178-182"
            , GantryAngle=178, ArcStopGantryAngle=182, ArcRotationDirection="CounterClockwise"
            , BeamQualityId="6"
            , IsocenterData=beamset.CreateDefaultIsocenterData(Position=isoCenter)   
            , CouchRotationAngle=0, CouchPitchAngle=0, CouchRollAngle=0, CollimatorAngle=20
        )
        beam.SetBolus(BolusName="")
        beam.SetDoseSpecificationPoint(Name='DSP')

        # Step 4 - Settings for dual arc beam
        plan.PlanOptimizations[0].OptimizationParameters.TreatmentSetupSettings[0].BeamSettings[0].ArcConversionPropertiesPerBeam.EditArcBasedBeamOptimizationSettings(CreateDualArcs=True
                                        , FinalGantrySpacing=2, MaxArcDeliveryTime=90, BurstGantrySpacing=None, MaxArcMU=None
                                    )
    except:
        traceback.print_exc()
        return False
    
    return True

########## Objective RELATED ##########

def getObjectivesFromPlan(plan, beamset):

    objective = plan.PlanOptimizations[plan.BeamSets.IndexOf(beamset)].Objective
    if objective is not None:
        objectives = objective.ConstituentFunctions
    else:
        objectives = []
    
    return objectives

def getObjectivesFromPath(plan, beamset, pathObjectives):
    
    templateManager = ObjectiveTemplateManager(plan, beamset)
    templateManager.parse_xml(str(pathObjectives))
    objectives = templateManager.objectives

    return objectives

def resetObjectives(plan, beamset):

    try:
        plan.PlanOptimizations[plan.BeamSets.IndexOf(beamset)].ClearConstituentFunctions() # plan.PlanOptimizations[0].Objective.ConstituentFunctions[0].DeleteFunction()
        constraintList = list(plan.PlanOptimizations[plan.BeamSets.IndexOf(beamset)].Constraints)
        for constraint in constraintList:
            constraint.DeleteFunction()
        
        return True
    
    except:
        traceback.print_exc()
        
    return False

def uploadObjectivesToRS(plan, beamset, objectivesFromPath):

    objectivesUploadStatus = False
    print ('\n - [uploadObjectivesToRS()] Uploading objectives to RS ...')

    try:
        # Step 1 - Reset objectives
        resetStatus = resetObjectives(plan, beamset)        
        if not resetStatus:
            print (' - [ERROR][uploadObjectivesToRS()] Could not reset objectives!')
            return objectivesUploadStatus
        
        # Step 2 - Upload objectives
        for index, objective in enumerate(objectivesFromPath):
            try:
                objective.apply()
                # print (f' - [uploadObjectivesToRS()] roi: {objective.roi_name}, fType: {objective.function_type}, weight: {objective.weight}')
            except:
                traceback.print_exc()
        
        objectivesUploadStatus = True

    except:
        traceback.print_exc()
        
    print (' - [uploadObjectivesToRS()] Uploaded objectives to RS: ', objectivesUploadStatus)
    return objectivesUploadStatus

########## Clinical Goals RELATED ##########

def getV95ForROI(case, beamset, roiName):
        
        v95 = None
        try:
            v95 = case.GetClinicalGoalValue(beamset, roiName, 'D95%')
        except:
            traceback.print_exc()
        
        return v95

########## Optimization RELATED ##########

def optimizePlan(planName, count=4, reset=True, pathIsoDoseXML=None):

    optimizeStatus = False
    objectiveValues = []
    try:
        # Step 1 - Get plan and beamset
        patient = rayStationSave()
        case    = patient.Cases[0]
        plan    = case.TreatmentPlans[planName]
        beamset = case.TreatmentPlans[planName].BeamSets[planName]
        beamSetIndex = plan.BeamSets.IndexOf(beamset)

        # Step 2 - Reset existing optimization
        if reset:
            plan.PlanOptimizations[beamSetIndex].ResetOptimization()
            print ('\n\n - [optimizePlan()][Patient={}][Plan={}] Optimization has been reset ...'.format(getPatientIdentifier(patient), planName))
        else:
            print ('\n\n - [optimizePlan()][Patient={}][Plan={}] Optimization has not been reset ...'.format(getPatientIdentifier(patient), planName))
            if plan.PlanOptimizations[beamSetIndex].ProgressOfOptimization is not None:
                objectiveValue = plan.PlanOptimizations[beamSetIndex].ProgressOfOptimization.ObjectiveValues[-1] # [TODO: range = (0.36,?)]
                print (' - [optimizePlan()][Patient={}][Plan={}] Previous objective {:.4f}'.format(getPatientIdentifier(patient), planName, objectiveValue))
        _ = rayStationSave()

        # Step 3 - Run optimization
        times = []
        for runID in range(count):
            t0 = time.time()
            print (' - [optimizePlan()][Patient={}][Plan={}] Running optimization step {}/{} ... '.format(getPatientIdentifier(patient), planName, runID+1, count))
            plan.PlanOptimizations[beamSetIndex].RunOptimization()   
            times.append(time.time() - t0)
            try:
                objectiveValue = plan.PlanOptimizations[beamSetIndex].ProgressOfOptimization.ObjectiveValues[-1] # [TODO: range = (0.36,?)]
                print (' --- [optimizePlan()] Optimization step {}/{} took {:.2f} seconds with mean objective: {:.4f}'.format(runID+1, count, times[-1], objectiveValue))
                objectiveValues.append(objectiveValue)
            except:
                traceback.print_exc()
        print (' - [optimizePlan()] Optimization took total {:.2f} seconds'.format(np.sum(times)))

        if runID == count-1:
            optimizeStatus = True
        
        if optimizeStatus:
            applyIsoDoseColors(pathIsoDoseXML, case, beamset)
        else:
            print (' - [optimizePlan()] Full optimization failed, not applying isodose colors')

    except:
        traceback.print_exc()
    
    return optimizeStatus, objectiveValues

########## RS ROI RELATED ##########

def checkOARDuplicateStatus(case, verbose=False):
    
    duplicateStatus = {}
    for oar in config.OARS:
        duplicateStatus[oar] = False
        if oar == config.KEYNAME_CAVITY_ORAL:
            oarDuplicate = config.KEYNAME_ORAL_CAVITY + config.KEY_AUTOCONTOUR_SUFFIX
        elif oar == config.KEYNAME_ESOPHAGUS_S:
            oarDuplicate = config.KEYNAME_ESOPHAGUS + config.KEY_AUTOCONTOUR_SUFFIX
        else:
            oarDuplicate = oar + config.KEY_AUTOCONTOUR_SUFFIX
        
        if config.RS_CHECK_PREFIX + oarDuplicate.translate(config.RS_STR_TRANSLATE_OBJ) in dir(case.PatientModel.StructureSets[0].RoiGeometries):
            duplicateStatus[oar] = True
            case.PatientModel.RegionsOfInterest[oarDuplicate].Color = config.OAR_DUPLICATE_COLOR_RGB_STRING # ["128, 0, 64", "Purple"]
        else:
            if verbose: print ('No duplicate found for', oar)
    
    return duplicateStatus

def checkROIExists(case, roiName):
    # [roi.Name for roi in case.PatientModel.RegionsOfInterest]
    exists = False

    roiNameCopy       = copy.deepcopy(roiName).replace('-(', '-') # 'Parotid_L-(CTV_DL1+3mm)'
    roiNameRSStandard = config.RS_CHECK_PREFIX + roiNameCopy.translate(config.RS_STR_TRANSLATE_OBJ)
    roiNameRSStandard = roiNameRSStandard.replace('___', '_')

    if roiNameRSStandard in dir(case.PatientModel.StructureSets[0].RoiGeometries):
        exists = True
    
    return exists

def getRoiVolume(case, roiName):

    volume = -1

    try:
        if checkROIExists(case, roiName):
            roiObj = case.PatientModel.StructureSets[0].RoiGeometries[roiName]
            if roiObj.HasContours():
                volume = roiObj.GetRoiVolume()
        else:
            print (f'\n - [getRoiVolume()] {roiName} does not exist in case: {case.CaseName}!')

    except Exception as e:
        print (e)
    
    return volume

# important function 1 (for project)
def doROIAlgebraForAutoContours(case):

    try:
        patientObj = connect.get_current(config.KEYNAME_PATIENT)
        patientID = patientObj.PatientID
    except:
        patientID = -1
        traceback.print_exc()

    def getMarginSettings(distInCm):
        return { 'Type': "Expand", 'Superior': distInCm, 'Inferior': distInCm, 'Anterior': distInCm, 'Posterior': distInCm, 'Right': distInCm, 'Left': distInCm }
    
    try:
        roiType = "Control"

        # For each ROI, create a new ROI with a 0.3cm margin
        if 1:
            params = [
                {
                    'roiNameNew': config.KEYNAME_BRAINSTEM + config.KEY_3MM_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX #'Brainstem+3 (1)'
                    , 'expARois': [config.KEYNAME_BRAINSTEM + config.KEY_AUTOCONTOUR_SUFFIX] # 'Brainstem (1)'
                    , 'expAMarginSettings': getMarginSettings(0.3)
                    , 'expBRois': []
                    , 'expBMarginSettings': getMarginSettings(0)
                }
                ,{
                    'roiNameNew': config.KEYNAME_SPINALCORD + config.KEY_3MM_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX #'SpinalCord+3 (1)'
                    , 'expARois': [config.KEYNAME_SPINALCORD + config.KEY_AUTOCONTOUR_SUFFIX] #'SpinalCord (1)'
                    , 'expAMarginSettings': getMarginSettings(0.3)
                    , 'expBRois': []
                    , 'expBMarginSettings': getMarginSettings(0)
                }
            ]

            for param in params:
                try:
                    roiName = param['expARois'][0]
                    roiNameNew = param['roiNameNew']
                    if checkROIExists(case, roiName) and not checkROIExists(case, roiNameNew):
                        case.PatientModel.CreateRoi(Name=roiNameNew, Color="128, 0, 64", Type=roiType, TissueName=None, RbeCellTypeName=None, RoiMaterial=None)
                        case.PatientModel.RegionsOfInterest[roiNameNew].SetAlgebraExpression(
                            ExpressionA={ 'Operation': "Union", 'SourceRoiNames': param['expARois'], 'MarginSettings': param['expAMarginSettings']}
                            , ExpressionB={ 'Operation': "Union", 'SourceRoiNames': param['expBRois'], 'MarginSettings': param['expBMarginSettings']}
                            , ResultOperation="None", ResultMarginSettings=getMarginSettings(0)
                        )
                        case.PatientModel.RegionsOfInterest[roiNameNew].Color = config.OAROBJ_DUPLICATE_COLOR_RGB_STRING 
                        case.PatientModel.RegionsOfInterest[roiNameNew].UpdateDerivedGeometry(Examination=case.Examinations['CT 1'], Algorithm="Auto")
                    else:
                        print (f' - [doROIAlgebraForAutoContours] {roiName}:{checkROIExists(case, roiName)} and {roiNameNew}:{checkROIExists(case, roiNameNew)}')
                
                except:
                    traceback.print_exc()

        # For each ROI, compare the ROI to PTV_DL1_DVH+0.5cm
        if 1:
            params = [
                {
                    'roiNameNew': config.KEYNAME_PAROTID_L + config.KEY_OBJ_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX #'Parotid_L_obj (1)'
                    , 'expARois': [config.KEYNAME_PAROTID_L + config.KEY_AUTOCONTOUR_SUFFIX] #'Parotid_L (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_PTV_DL1_DVH] #'PTV_DL1_DVH'
                    , 'expBMarginSettings': getMarginSettings(0.5)
                }
                ,{
                    'roiNameNew': config.KEYNAME_PAROTID_R + config.KEY_OBJ_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX #'Parotid_R_obj (1)'
                    , 'expARois': [config.KEYNAME_PAROTID_R + config.KEY_AUTOCONTOUR_SUFFIX] #'Parotid_R (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_PTV_DL1_DVH] #'PTV_DL1_DVH'
                    , 'expBMarginSettings': getMarginSettings(0.5)
                }
                ,{
                    'roiNameNew': config.KEYNAME_ORAL_CAVITY + config.KEY_OBJ_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX #'Oral_Cavity_obj (1)'
                    , 'expARois': [config.KEYNAME_ORAL_CAVITY + config.KEY_AUTOCONTOUR_SUFFIX] # 'Oral_Cavity (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_PTV_DL1_DVH] #'PTV_DL1_DVH'
                    , 'expBMarginSettings': getMarginSettings(0.5)
                }
                , {
                    'roiNameNew': config.KEYNAME_SUBMAND_L_OBJ + config.KEY_AUTOCONTOUR_SUFFIX # 'Submand_L_obj (1)'
                    , 'expARois': [config.KEYNAME_GLND_SUBMAND_L + config.KEY_AUTOCONTOUR_SUFFIX] # 'Glnd_Submand_L (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_PTV_DL1_DVH] #'PTV_DL1_DVH'
                    , 'expBMarginSettings': getMarginSettings(0.5)
                }
                , {
                    'roiNameNew': config.KEYNAME_SUBMAND_R_OBJ + config.KEY_AUTOCONTOUR_SUFFIX # 'Submand_R_obj (1)'
                    , 'expARois': [config.KEYNAME_GLND_SUBMAND_R + config.KEY_AUTOCONTOUR_SUFFIX] # 'Glnd_Submand_R (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_PTV_DL1_DVH] #'PTV_DL1_DVH'
                    , 'expBMarginSettings': getMarginSettings(0.5)
                }
                , {
                    'roiNameNew': config.KEYNAME_MANDIBLE_PTV + config.KEY_AUTOCONTOUR_SUFFIX # 'Bone_Mandible-PTV_DL1 (1)'
                    , 'expARois': [config.KEYNAME_MANDIBLE + config.KEY_AUTOCONTOUR_SUFFIX] # 'Glnd_Submand_R (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_PTV_DL1_DVH] #'PTV_DL1_DVH'
                    , 'expBMarginSettings': getMarginSettings(0.0)
                }
            ]

            # For params[-1] i.e. Bone_Mandible-PTV_DL1
            if 1:
                if patientID in ['HCAI-Dose-x5']:
                    print (' - [doROIAlgebraForAutoContours] For Bone_Mandible-PTV_DL1 doing getMarginSettings(0.5) for patient: ', patientID)
                    params[-1]['expBMarginSettings'] = getMarginSettings(0.5)
            
            for param in params:
                try:
                    roiNameNew = param['roiNameNew']
                    roiName    = param['expARois'][0]
                    roiNameRef = param['expBRois'][0]
                    if checkROIExists(case, roiName) and checkROIExists(case, roiNameRef) and not checkROIExists(case, roiNameNew):
                        case.PatientModel.CreateRoi(Name=roiNameNew, Color="Red", Type=roiType, TissueName=None, RbeCellTypeName=None, RoiMaterial=None)
                        case.PatientModel.RegionsOfInterest[roiNameNew].SetAlgebraExpression(
                            ExpressionA={ 'Operation': "Union", 'SourceRoiNames': param['expARois'], 'MarginSettings': param['expAMarginSettings'] }
                            , ExpressionB={ 'Operation': "Union", 'SourceRoiNames': param['expBRois'], 'MarginSettings': param['expBMarginSettings']}
                            , ResultOperation="Subtraction", ResultMarginSettings=getMarginSettings(0)
                        )
                        case.PatientModel.RegionsOfInterest[roiNameNew].Color = config.OAROBJ_DUPLICATE_COLOR_RGB_STRING 
                        case.PatientModel.RegionsOfInterest[roiNameNew].UpdateDerivedGeometry(Examination=case.Examinations['CT 1'], Algorithm="Auto")
                    else:
                        print (f' - [doROIAlgebraForAutoContours] {roiName}:{checkROIExists(case, roiName)} and {roiNameNew}:{checkROIExists(case, roiNameNew)} and {roiNameRef}:{checkROIExists(case, roiNameRef)}')

                except:
                    traceback.print_exc()

        # For swallowing muscles
        if 1:
            try:
                roiNameNew = config.KEYNAME_SWAL_COMP + config.KEY_AUTOCONTOUR_SUFFIX # 'Swal_Comp (1)'
                if not checkROIExists(case, roiNameNew):
                    roiNew = case.PatientModel.CreateRoi(Name=roiNameNew, Color=config.OAR_DUPLICATE_COLOR_RGB_STRING, Type=roiType, TissueName=None, RbeCellTypeName=None, RoiMaterial=None)
                    roiName1 = config.KEYNAME_MUSC_CONSTRICT_I
                    roiName2 = config.KEYNAME_MUSC_CONSTRICT_M
                    roiName3 = config.KEYNAME_MUSC_CONSTRICT_S
                    roiName4 = config.KEYNAME_CRICO
                    roiName5 = config.KEYNAME_LARYNX_SG + config.KEY_AUTOCONTOUR_SUFFIX # 'Larynx_SG (1)'
                    roiName6 = config.KEYNAME_GLOTTIC_AREA

                    sourceRoiNamesPotential = [roiName1, roiName2, roiName3, roiName4, roiName5, roiName6]
                    sourceRoiNames = []
                    for sourceRoiName in sourceRoiNamesPotential:
                        if checkROIExists(case, sourceRoiName):
                            sourceRoiNames.append(sourceRoiName)
                        else:
                            print (f' - [doROIAlgebraForAutoContours] {sourceRoiName} not available!')

                    if len(sourceRoiNames):
                        roiNew.SetAlgebraExpression(
                            ExpressionA={ 'Operation': "Union", 'SourceRoiNames': sourceRoiNames, 'MarginSettings': getMarginSettings(0)}
                            , ExpressionB={ 'Operation': "Union", 'SourceRoiNames': [], 'MarginSettings': getMarginSettings(0)}
                            , ResultOperation="None", ResultMarginSettings=getMarginSettings(0)
                        )
                        roiNew.UpdateDerivedGeometry(Examination=case.Examinations['CT 1'], Algorithm="Auto")
                        case.PatientModel.RegionsOfInterest[roiNameNew].Color = config.OAROBJ_DUPLICATE_COLOR_RGB_STRING 
                        case.PatientModel.RegionsOfInterest[roiNameNew].UpdateDerivedGeometry(Examination=case.Examinations['CT 1'], Algorithm="Auto")

                    else:
                        print (f' - [doROIAlgebraForAutoContours] {roiName1}:{checkROIExists(case, roiName1)} and {roiName2}:{checkROIExists(case, roiName2)} and {roiName3}:{checkROIExists(case, roiName3)} and {roiName4}:{checkROIExists(case, roiName4)} and {roiName5}:{checkROIExists(case, roiName5)} and {roiName6}:{checkROIExists(case, roiName6)}')
                else:
                    print (f' - [doROIAlgebraForAutoContours] {roiNameNew}:{checkROIExists(case, roiNameNew)}')

                roiNameNewObj = config.KEYNAME_SWAL_OBJ + config.KEY_AUTOCONTOUR_SUFFIX # 'Swal_obj (1)'
                if checkROIExists(case, roiNameNew) and not checkROIExists(case, roiNameNewObj):
                    case.PatientModel.CreateRoi(Name=roiNameNewObj, Color=config.OAR_DUPLICATE_COLOR_RGB_STRING, Type=roiType, TissueName=None, RbeCellTypeName=None, RoiMaterial=None)

                    if checkROIExists(case, config.KEYNAME_PTV_DL2_DVH):
                        case.PatientModel.RegionsOfInterest[roiNameNewObj].SetAlgebraExpression(
                            ExpressionA={ 'Operation': "Union", 'SourceRoiNames': [config.KEYNAME_SWAL_COMP + config.KEY_AUTOCONTOUR_SUFFIX], 'MarginSettings': getMarginSettings(0)}
                            , ExpressionB={ 'Operation': "Union", 'SourceRoiNames': [config.KEYNAME_PTV_DL1_DVH, config.KEYNAME_PTV_DL2_DVH], 'MarginSettings': getMarginSettings(0.5)}
                            , ResultOperation="Subtraction", ResultMarginSettings=getMarginSettings(0)
                        )
                    else:
                        case.PatientModel.RegionsOfInterest[roiNameNewObj].SetAlgebraExpression(
                            ExpressionA={ 'Operation': "Union", 'SourceRoiNames': [config.KEYNAME_SWAL_COMP + config.KEY_AUTOCONTOUR_SUFFIX], 'MarginSettings': getMarginSettings(0)}
                            , ExpressionB={ 'Operation': "Union", 'SourceRoiNames': [config.KEYNAME_PTV_DL1_DVH], 'MarginSettings': getMarginSettings(0.5)}
                            , ResultOperation="Subtraction", ResultMarginSettings=getMarginSettings(0)
                        )
                    case.PatientModel.RegionsOfInterest[roiNameNewObj].Color = config.OAROBJ_DUPLICATE_COLOR_RGB_STRING 
                    case.PatientModel.RegionsOfInterest[roiNameNewObj].UpdateDerivedGeometry(Examination=case.Examinations['CT 1'], Algorithm="Auto")
                else:
                    print (f' - [doROIAlgebraForAutoContours] {roiNameNew}:{checkROIExists(case, roiNameNew)} and {roiNameNewObj}:{checkROIExists(case, roiNameNewObj)}')
            except:
                traceback.print_exc()

    
    except:
        traceback.print_exc()

def doROIAlgebraForProtonAutoContours(case):
    """
    1. opt_CTV_L, opt_CTV_R (depends on parotid L/R)
        - Box_L = min_z(Parotid_L) --> min_z(CTV_DL1)
        - Box_R = min_z(Parotid_R) --> min_z(CTV_DL1)
        * opt_CTV_L = CTV_DL1 - (Box_R)
        * opt_CTV_R = CTV_DL1 - (Box_L)
    2. _obj like rois
        - Oral_Cavity-(CTV_DL1+0.3cm)
    
    NOTE: This is NOT a force operation (i.e. does not do anything if the roi already exists) 
    """

    def getMarginSettings(distInCm):
        return { 'Type': "Expand", 'Superior': distInCm, 'Inferior': distInCm, 'Anterior': distInCm, 'Posterior': distInCm, 'Right': distInCm, 'Left': distInCm }
    
    try:
        roiType = "Control"
        examinationName = case.Examinations[0].Name

        # For a ROI, compare the ROI to CTV_DL1+0.3cm
        if 1:
            params = [
                {
                    'roiNameNew': config.KEYNAME_PAROTID_L + config.KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX # 'Parotid_L-(CTV_DL1+3mm) (1)'
                    , 'expARois': [config.KEYNAME_PAROTID_L + config.KEY_AUTOCONTOUR_SUFFIX] #'Parotid_L (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_ROI_CTV_DL1] # 'CTV_DL1'
                    , 'expBMarginSettings': getMarginSettings(0.3)
                }
                , {
                    'roiNameNew': config.KEYNAME_PAROTID_R + config.KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX #'Parotid_R-(CTV_DL1+3mm) (1)'
                    , 'expARois': [config.KEYNAME_PAROTID_R + config.KEY_AUTOCONTOUR_SUFFIX] #'Parotid_R (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_ROI_CTV_DL1] # 'CTV_DL1'
                    , 'expBMarginSettings': getMarginSettings(0.3)
                }
                , {
                    'roiNameNew': config.KEYNAME_ORAL_CAVITY + config.KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX # Oral_Cavity-(CTV_DL1+3mm) (1)
                    , 'expARois': [config.KEYNAME_ORAL_CAVITY + config.KEY_AUTOCONTOUR_SUFFIX] # 'Oral_Cavity (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_ROI_CTV_DL1] # 'CTV_DL1'
                    , 'expBMarginSettings': getMarginSettings(0.3)
                }
                , {
                    'roiNameNew': config.KEYNAME_GLND_SUBMAND_L + config.KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX # 'Submand_L_obj-(CTV_DL1+3mm) (1)'
                    , 'expARois': [config.KEYNAME_GLND_SUBMAND_L + config.KEY_AUTOCONTOUR_SUFFIX] # 'Glnd_Submand_L (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_ROI_CTV_DL1] # 'CTV_DL1'
                    , 'expBMarginSettings': getMarginSettings(0.3)
                }
                , {
                    'roiNameNew': config.KEYNAME_GLND_SUBMAND_R + config.KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX # 'Submand_R_obj-(CTV_DL1+3mm) (1)'
                    , 'expARois': [config.KEYNAME_GLND_SUBMAND_R + config.KEY_AUTOCONTOUR_SUFFIX] # 'Glnd_Submand_R (1)'
                    , 'expAMarginSettings': getMarginSettings(0)
                    , 'expBRois': [config.KEYNAME_ROI_CTV_DL1] # 'CTV_DL1'
                    , 'expBMarginSettings': getMarginSettings(0.3)
                }
            ]

            for param in params:
                try:
                    roiNameNew = param['roiNameNew']
                    roiName    = param['expARois'][0]
                    roiNameRef = param['expBRois'][0]
                    if checkROIExists(case, roiName) and checkROIExists(case, roiNameRef) and getRoiVolume(case, roiNameNew) <= 0.0:
                        if not checkROIExists(case, roiNameNew):
                            case.PatientModel.CreateRoi(Name=roiNameNew, Color="Red", Type=roiType, TissueName=None, RbeCellTypeName=None, RoiMaterial=None)
                        case.PatientModel.RegionsOfInterest[roiNameNew].SetAlgebraExpression(
                            ExpressionA={ 'Operation': "Union", 'SourceRoiNames': param['expARois'], 'MarginSettings': param['expAMarginSettings'] }
                            , ExpressionB={ 'Operation': "Union", 'SourceRoiNames': param['expBRois'], 'MarginSettings': param['expBMarginSettings']}
                            , ResultOperation="Subtraction", ResultMarginSettings=getMarginSettings(0)
                        )
                        case.PatientModel.RegionsOfInterest[roiNameNew].Color = config.OAROBJ_DUPLICATE_COLOR_RGB_STRING 
                        case.PatientModel.RegionsOfInterest[roiNameNew].UpdateDerivedGeometry(Examination=case.Examinations[0], Algorithm="Auto")
                        assert getRoiVolume(case, roiNameNew) > 0.0, f' - [doROIAlgebraForProtonAutoContours] No volume for {roiNameNew}:{checkROIExists(case, roiNameNew)}'
                    else:
                        print (f'  -- [doROIAlgebraForProtonAutoContours] {roiName}:{checkROIExists(case, roiName)} and {roiNameRef}:{checkROIExists(case, roiNameRef)} and {roiNameNew}:{checkROIExists(case, roiNameNew), getRoiVolume(case, roiNameNew)}')

                except:
                    print (f' - [doROIAlgebraForProtonAutoContours] {roiName}:{checkROIExists(case, roiName)} and {roiNameNew}:{checkROIExists(case, roiNameNew)} and {roiNameRef}:{checkROIExists(case, roiNameRef)}')
                    traceback.print_exc()

        # For Mid_structures
        # (esophagus (1) + trachea + larynx_sg (1) + glottic_area) - CTV_DL1+3mm = Mid_Structures-(CTV_DL1+3mm) (1)
        if 1:
            try:
                roiNameNew = config.KEYNAME_ROI_MID_STRUCTURES + config.KEYNAME_ROI_CTV_DL1_3MM_AS_SUFFIX + config.KEY_AUTOCONTOUR_SUFFIX # 'Mid_Structures-(CTV_DL1+3mm) (1)'
                if getRoiVolume(case, roiNameNew) <= 0.0:
                    if not checkROIExists(case, roiNameNew):
                        case.PatientModel.CreateRoi(Name=roiNameNew, Color=config.OAR_DUPLICATE_COLOR_RGB_STRING, Type=roiType, TissueName=None, RbeCellTypeName=None, RoiMaterial=None)
                    roiName1 = config.KEYNAME_ESOPHAGUS + config.KEY_AUTOCONTOUR_SUFFIX # 'Esophagus (1)'
                    roiName2 = config.KEYNAME_TRACHEA 
                    roiName3 = config.KEYNAME_LARYNX_SG + config.KEY_AUTOCONTOUR_SUFFIX # 'Larynx_SG (1)'
                    roiName4 = config.KEYNAME_GLOTTIC_AREA

                    sourceRoiNamesPotential = [roiName1, roiName2, roiName3, roiName4]
                    sourceRoiNames = []
                    for sourceRoiName in sourceRoiNamesPotential:
                        if checkROIExists(case, sourceRoiName):
                            sourceRoiNames.append(sourceRoiName)
                        else:
                            print (f' - [doROIAlgebraForProtonAutoContours] {sourceRoiName} not available!')

                    if len(sourceRoiNames):
                        case.PatientModel.RegionsOfInterest[roiNameNew].SetAlgebraExpression(
                                ExpressionA={ 'Operation': "Union", 'SourceRoiNames': sourceRoiNames, 'MarginSettings': getMarginSettings(0)}
                                , ExpressionB={ 'Operation': "Union", 'SourceRoiNames': [config.KEYNAME_ROI_CTV_DL1], 'MarginSettings': getMarginSettings(0.3)}
                                , ResultOperation="Subtraction", ResultMarginSettings=getMarginSettings(0)
                            )
                        case.PatientModel.RegionsOfInterest[roiNameNew].UpdateDerivedGeometry(Examination=case.Examinations[0], Algorithm="Auto")
                        case.PatientModel.RegionsOfInterest[roiNameNew].Color = config.OAROBJ_DUPLICATE_COLOR_RGB_STRING 
                        case.PatientModel.RegionsOfInterest[roiNameNew].UpdateDerivedGeometry(Examination=case.Examinations[0], Algorithm="Auto")
                        assert getRoiVolume(case, roiNameNew) > 0.0, f' - [doROIAlgebraForProtonAutoContours] No volume for {roiNameNew}:{checkROIExists(case, roiNameNew)}'

                    else:
                        print (f' - [doROIAlgebraForAutoContours] {roiName1}:{checkROIExists(case, roiName1)} and {roiName2}:{checkROIExists(case, roiName2)} and {roiName3}:{checkROIExists(case, roiName3)} and {roiName4}:{checkROIExists(case, roiName4)}')
                else:
                    print (f' - [doROIAlgebraForProtonAutoContours] {roiNameNew}:{checkROIExists(case, roiNameNew)}, volume={getRoiVolume(case, roiNameNew)}')
            
            except:
                print (f' - [doROIAlgebraForProtonAutoContours] {roiName1}:{checkROIExists(case, roiName1)} and {roiName2}:{checkROIExists(case, roiName2)} and {roiName3}:{checkROIExists(case, roiName3)} and {roiName4}:{checkROIExists(case, roiName4)}')
                traceback.print_exc()

    except:
        traceback.print_exc()

def getRoiRootFromXML(pathXMLObj):

    try:
        
        tree = ElementTree.parse(str(pathXMLObj))
        root = tree.getroot()
        if root is None:
            print(f' - [updateKNOXMLForAutoContouring()] {config.KEY_XML_ROOT_NAME} does not exist')
            return tree, None

        subRoot = root.find(config.KEY_XML_SUBROOT_NAME)
        if subRoot is None:
            print(f' - [updateKNOXMLForAutoContouring()] {config.KEY_XML_SUBROOT_NAME} does not exist')
            return tree, None

        return tree, subRoot

    except:
        traceback.print_exc()
    
    return tree, None

# important function 2 (for project)
def updateKNOXMLForAutoContours(pathInitialObj, pathPatientObj, potentialRoisToRenameInAuto):

    # Step 1 - check if OG KNO.xml exists
    pathCSAutoObj = None
    if not Path(pathInitialObj).exists():
        print (f' - [updateKNOXMLForAutoContouring()] Initial objectives XML file not found at {pathInitialObj}')
        return None
    if not Path(pathPatientObj).exists():
        print (f' - [updateKNOXMLForAutoContouring()] Patient-Solution XML file not found at {pathPatientObj}')
        return None
    
    # Step 2 - Create a list of ROI names to replace
    finalRoisToRenameInAuto = []
    print (f'\n - [updateKNOXMLForAutoContouring()] POTENTIALLY renaming these ROIs in auto xml file: {potentialRoisToRenameInAuto}')

    # Step 3 - Read pathKNOPatientObj and check if it contains roi from potentialRoisToRenameInAuto
    _, patientObjRoiRoot = getRoiRootFromXML(pathPatientObj)
    for roiEl in patientObjRoiRoot.getchildren():
        try:
            if roiEl.tag == config.KEY_XML_ROI:
                if roiEl.get(config.KEY_XML_NAME) in potentialRoisToRenameInAuto:
                    finalRoisToRenameInAuto.append(roiEl.get(config.KEY_XML_NAME))
        except:
            traceback.print_exc()
    print (f' - [updateKNOXMLForAutoContouring()] FINALLY renaming these ROIs in auto xml file: {finalRoisToRenameInAuto}')

    # Step 4 - Rename rois from finalRoisToRenameInAuto in pathKNOCSObj
    initialObjTree, initialObjRoiRoot = getRoiRootFromXML(pathInitialObj)
    for roiEle in initialObjRoiRoot.getchildren():
        try:
            
            if roiEle.tag == config.KEY_XML_ROI:
                roiElName = roiEle.get(config.KEY_XML_NAME)
                # print (f' - [updateKNOXMLForAutoContouring()] {roiElName} in finalRoisToRenameInAuto = { roiEle.get(config.KEY_XML_NAME) in finalRoisToRenameInAuto} ')
                if roiElName in finalRoisToRenameInAuto:
                    roiEle.set(config.KEY_XML_NAME, roiElName + config.KEY_AUTOCONTOUR_SUFFIX)
        except:
            traceback.print_exc()

    # Step 4 - Save the XML file
    pathAutoObj = Path(pathInitialObj).parent.absolute().joinpath(Path(pathInitialObj).stem + config.KEY_AUTOCONTOUR_SUFFIX + Path(pathInitialObj).suffix)
    initialObjTree.write(str(pathAutoObj))

    return pathAutoObj

##########################################################################################
#                                     FRANKS CODE                                        #
##########################################################################################

# https://git.lumc.nl/fjwmdankers/rs-batch-dvh-parameters/functions_dvh_parameters.py
def read_dvhparamlist_csv_or_txt(csv_path, csv_delimiter):
    # reads roi/dvh parameters from a csv/txt file
    # 
    # Input data in csv should be formatted such that:
    #   [roiname1];[dvh param 1];[dvh param 2][dvh param 3]
    #   [roiname2];[dvh param 1]
    #
    # Notes on the DVH param notation:
    #   "Volume" must have no output unit (cc is assumed)
    #   "Dmax (Gy)", "Dmin (Gy)", "Dmean (Gy)" must have an output unit (cGy/Gy)
    #   relative doses (input or output) should also specify the doselevel (e.g., "D2cc (%DL1)", "V95%DL1 (cc)", "CI95%DL1 (RTOG)")
    #
    # Output of this functions looks like this for example:
    #    rois_and_dvhparams = {
    #           "CTV_DL2_DVH": ["V99%DL2 (%)","V99%DL2 (cc)","V20Gy (cc)","V2000cGy (%)","Volume"],
    #           "CTV_DL1_DVH": ["D2cc (%DL2)","D90% (%DL1)","D90% (cGy)","D90% (Gy)", "Dmax (cGy)", "Dmean (Gy)", "Dmin (%DL1)"],
    #           "PTV_DL1_DVH": ["Volume", "V95%DL1 (%)", "CI95%DL1 (Riet)", "CI95%DL1 (RTOG)", "CI95%DL1 (RS)", "HI95%"]
    #   }

    # open file
    csv = open(csv_path, 'r')

    # loop as long as there are new lines in the csv
    count = 0
    rois_and_dvhparams = {}
    while True:
        # read line, check if we're done, and increase counter
        line = csv.readline()
        if not line:
            break
        count += 1

        # some csv have encodings that lead to entry characters at the first line
        line = line.replace('ï»¿','')
        line = line.replace('\n','')

        # process the line
        rowlist = line.split(csv_delimiter)
        
        # remove empty entries (e.g., due to trailing delimiters)
        rowlist = [x for x in rowlist if str(x) != '']
        
        # save to dictionary
        if len(rowlist) > 0: # only add if the row is not empty, e.g., there might be empty rows in excel
            if rowlist[0] not in rois_and_dvhparams:
                rois_and_dvhparams[rowlist[0]] = rowlist[1:]
            else:
                rois_and_dvhparams[rowlist[0]] = rois_and_dvhparams[rowlist[0]] + rowlist[1:]

    return(rois_and_dvhparams)

def process_dvhparam(plan, doselevels, dose, roi, dvhparam_label, verbose=False):
    
    # dvhparam not recognized, return empty string (leads to empty entry in csv)
    dvhparam_quantity, dvhparam_unit, dvhparam_input_value, dvhparam_input_unit, dvhparam_untangle_status = untangle_dvhparam_string(dvhparam_label, verbose=verbose)

    if dvhparam_untangle_status:
        dvhparam_value = calc_dvhparam(plan, doselevels, dose, roi, dvhparam_label, dvhparam_quantity, dvhparam_unit, dvhparam_input_value, dvhparam_input_unit, verbose)
        if verbose: print(f"\t\t\t{round(dvhparam_value,5)}")
    else:
        # dvhparam label untangling failed
        dvhparam_value = ''

    return(dvhparam_value)

def calc_dvhparam(plan, doselevels, dose, roi, dvhparam_label, dvhparam_quantity, dvhparam_unit, dvhparam_input_value, dvhparam_input_unit, verbose=False):
    # TODO add more parameter processing options

    try:
        dvhparam_value = ''
        exam = connect.get_current("Examination")
        case = connect.get_current("Case")

        # derive dvh parameter value based on the desired combinations of input/output units etc
        if dvhparam_quantity == 'V':
            if dvhparam_label == 'Volume':
                # volume = dose.GetDoseGridRoi(RoiName=roi).RoiVolumeDistribution.TotalVolume             # based on plan/fraction dose (this is what dose statistics table also shows, and roi properties is close to this)
                # beamset = plan.BeamSets[0]
                volume = case.PatientModel.StructureSets[exam.Name].RoiGeometries[roi].GetRoiVolume()  # based on structureset (this value is typically slightly larger than what RS shows in several windows, e.g. roi properties, dose statistics)
                # volume3 = beamset.GetStructureSet().RoiGeometries[roi].GetRoiVolume()                   # same result as above
                # print("----------------------")
                # print(f"volume based on dosegrid:  {volume}")
                # print(f"volume based on exam:  {volume2}")
                # print(f"volume based on beamset:  {volume3}")
                dvhparam_value = volume
            elif dvhparam_input_unit in ['Gy', 'cGy', '%DL1', '%DL2', '%DL3', '%DL4', '%DL5', '%DL6']:
                # process dvh input value using the dvh input unit
                if dvhparam_input_unit == 'Gy':
                    dvhparam_input_value = dvhparam_input_value*100
                elif dvhparam_input_unit == 'cGy':
                    pass # we need cGy so nothing has to be done
                elif '%DL' in dvhparam_input_unit:
                    doselevel_ref_idx = int(dvhparam_input_unit[-1])-1
                    dvhparam_input_value = dvhparam_input_value / 100 * doselevels[doselevel_ref_idx]

                # get dvh parameter output value
                dvhparam_value_rel = dose.GetRelativeVolumeAtDoseValues(RoiName=roi, DoseValues=[dvhparam_input_value])[0] # fraction of total volume

                # prepare for desired output unit
                if dvhparam_unit == '%':
                    dvhparam_value = dvhparam_value_rel * 100
                elif dvhparam_unit == 'cc':
                    volume = dose.GetDoseGridRoi(RoiName=roi).RoiVolumeDistribution.TotalVolume                     # based on plan/fraction dose (this is what dose statistics table also shows, and roi properties is close to this)
                    #volume = case.PatientModel.StructureSets[exam.Name].RoiGeometries['Lung_R'].GetRoiVolume()     # based on structureset (this value is typically larger than what RS shows in several windows, e.g. roi properties, dose statistics)
                    dvhparam_value = dvhparam_value_rel * volume
                else:
                    dvhparam_value = ''
                    print(f"\t\tcant process this dvhparam_unit yet ({dvhparam_unit})")
            else:
                print(f"\t\tcant process this dvhparam_input_unit yet ({dvhparam_input_unit})")
        elif dvhparam_quantity == 'D':
            if 'Dmin' in dvhparam_label or 'Dmax' in dvhparam_label or 'Dmean' in dvhparam_label:
                if 'Dmin' in dvhparam_label:
                    dvhparam_value = dose.GetDoseStatistic(RoiName=roi, DoseType='Min')                             # equivalent to:    dose.GetDoseAtRelativeVolumes(RoiName=roiname, RelativeVolumes=[1])[0]
                elif 'Dmax' in dvhparam_label:
                    dvhparam_value = dose.GetDoseStatistic(RoiName=roi, DoseType='Max')                             # equivalent to:    dose.GetDoseAtRelativeVolumes(RoiName=roiname, RelativeVolumes=[0.00000001/volume])[0]
                elif 'Dmean' in dvhparam_label:
                    dvhparam_value = dose.GetDoseStatistic(RoiName=roi, DoseType='Average')
            elif dvhparam_input_unit in ['cc', '%']:
                # process dvh input value using the dvh input unit
                if dvhparam_input_unit == 'cc':
                    volume = dose.GetDoseGridRoi(RoiName=roi).RoiVolumeDistribution.TotalVolume                     # based on plan/fraction dose (this is what dose statistics table also shows, and roi properties is close to this)
                    #volume = case.PatientModel.StructureSets[exam.Name].RoiGeometries['Lung_R'].GetRoiVolume()     # based on structureset (this value is typically larger than what RS shows in several windows, e.g. roi properties
                    dvhparam_input_value = dvhparam_input_value/volume
                else:
                    dvhparam_input_value = dvhparam_input_value / 100 # needs to be relative between 0 and 1

                # get dvh parameter output value
                dvhparam_value = dose.GetDoseAtRelativeVolumes(RoiName=roi, RelativeVolumes=[dvhparam_input_value])[0]
                
            # prepare for desired output
            if '%DL' in dvhparam_unit:
                doselevel_ref_idx = int(dvhparam_unit[-1])-1
                dvhparam_value = dvhparam_value / doselevels[doselevel_ref_idx] * 100
            elif dvhparam_unit == 'Gy':
                dvhparam_value = dvhparam_value / 100
            elif dvhparam_unit == 'cGy':
                pass # it's already in cGy so nothing has to be done
            else:
                dvhparam_value = ''
                print(f"\t\tcant process this dvhparam_unit yet ({dvhparam_unit})")
        elif dvhparam_quantity == 'CI':
            # Raystation uses Loman/Scheib definition:
            #   CI_raystation = TVri / Vri
            #   CI_rtog = Vri / TV
            #   CI_riet = TVri2 / (TV * Vri)
            #
            # with:
            #   Vri = volume of reference isodose
            #   TV = target volume (most often the PTV volume)
            #   TVri = PTV volume covered by the refence isodose
            
            # find and set the external ROI
            # note: the external is used for a bugfix (forgot which one), if it doesnt find the external then we can't calculate the CI at the moment
            roi_external = ''
            roi_external_options = ['External', 'Body'] #hard-coded list for now, External=Pinnacle era, Body=RayStation era
            rois_case = [roi.Name for roi in case.PatientModel.RegionsOfInterest]
            for roi_external_option in roi_external_options:
                if roi_external_option in rois_case:
                    roi_external = roi_external_option
                    print(f"\t\t\tbody contour used for calculation: {roi_external}")

            # process dvh input value to cGy using the dvh input unit
            if dvhparam_input_unit == 'Gy':
                dvhparam_input_value = dvhparam_input_value*100
            elif dvhparam_input_unit == 'cGy':
                pass # we need cGy so nothing has to be done
            elif '%DL' in dvhparam_input_unit:
                doselevel_ref_idx = int(dvhparam_input_unit[-1])-1
                dvhparam_input_value = dvhparam_input_value / 100 * doselevels[doselevel_ref_idx]
            
            if roi_external:
                # calculate the parameters of the CI formulas (Vri, TV, TVri)
                external_volume = dose.GetDoseGridRoi(RoiName=roi_external).RoiVolumeDistribution.TotalVolume
                external_volume_covered_by_ref_isodose_relative = dose.GetRelativeVolumeAtDoseValues(RoiName=roi_external, DoseValues=[dvhparam_input_value])[0]
                volume_ref_isodose = external_volume_covered_by_ref_isodose_relative * external_volume          # using external volume to get the volume ref isodose (was een bugfix ergens voor, voor als contour buiten external zit)
            else:
                print(f"\t\t\texternal/body ROI not detected, calulated CI without, ~1 prct different (checked for: {', '.join(roi_external_options)})")
                roi_ref_temp_name = 'isodose_temp'
                roi_ref_temp = case.PatientModel.CreateRoi(Name=roi_ref_temp_name, Color='red', Type='Control')
                roi_ref_temp.CreateRoiGeometryFromDose(DoseDistribution=dose, ThresholdLevel=int(dvhparam_input_value))
                volume_ref_isodose = case.PatientModel.StructureSets[exam.Name].RoiGeometries[roi_ref_temp_name].GetRoiVolume()    # based on structureset (this value is typically slightly larger than what RS shows in several windows, e.g. roi properties, dose statistics), we use this here so we don't need to update the dose statistics (which is not possible if the machine is deprecated, and potentially if the plan is approved)
                case.PatientModel.RegionsOfInterest[roi_ref_temp_name].DeleteRoi()

            target_volume = case.PatientModel.StructureSets[exam.Name].RoiGeometries[roi].GetRoiVolume()    # based on structureset (this value is typically slightly larger than what RS shows in several windows, e.g. roi properties, dose statistics)
            # target_volume = dose.GetDoseGridRoi(RoiName=roi).RoiVolumeDistribution.TotalVolume            # based on plan/fraction dose (this is what dose statistics table also shows, and roi properties is close to this)
            target_volume_covered_by_ref_isodose_relative = dose.GetRelativeVolumeAtDoseValues(RoiName=roi, DoseValues=[dvhparam_input_value])[0]

            # calculate the CI itself
            if volume_ref_isodose and target_volume:
                if dvhparam_unit == "RTOG":
                    dvhparam_value = (volume_ref_isodose / target_volume)
                elif dvhparam_unit == "Riet":
                    target_volume_covered_by_ref_isodose = target_volume_covered_by_ref_isodose_relative * target_volume
                    dvhparam_value = ((target_volume_covered_by_ref_isodose ** 2) / (target_volume * volume_ref_isodose))
                elif dvhparam_unit == "RS":
                    target_volume_covered_by_ref_isodose = target_volume_covered_by_ref_isodose_relative * target_volume
                    dvhparam_value = target_volume_covered_by_ref_isodose / volume_ref_isodose
            else:
                print(f"\t\t\tcant calculate CI, either volume_ref_isodose ({volume_ref_isodose}) or roi_volume ({target_volume}) missing")

        elif dvhparam_quantity == 'HI':
            # Di%/D(1-i)%
            dvhparam_input_value = dvhparam_input_value / 100
            dose_at_parameter_volume = dose.GetDoseAtRelativeVolumes(RoiName=roi, RelativeVolumes=[dvhparam_input_value])[0]
            dose_at_remaining_volume = dose.GetDoseAtRelativeVolumes(RoiName=roi, RelativeVolumes=[1.0 - dvhparam_input_value])[0]

            dvhparam_value = dose_at_parameter_volume / dose_at_remaining_volume

        else:
            # can't process this dvhparam (yet)
            print(f"\t\t - [ERROR] cant process this dvhparam_quantity yet ({dvhparam_quantity})")
    except:
        print(f"\t\t - [ERROR] cant process this dvhparam_label ({dvhparam_label})")
        traceback.print_exc()
        

    return(dvhparam_value)

def untangle_dvhparam_string(dvhparam_label, verbose=False):
    output_quantity = ''
    output_unit = ''
    input_value = ''
    input_unit = ''
    untangle_status = ''

    # untangle the dvh parameter string using regexp
    if verbose: print(f"\t\t{dvhparam_label}")
    if dvhparam_label == 'Volume':
        output_quantity = 'V'
        output_unit = 'cc'
        untangle_status = 1
    elif 'Dmin' in dvhparam_label or 'Dmax' in dvhparam_label or 'Dmean' in dvhparam_label:
        output_quantity = 'D'
        output_unit = re.search('\((.*)\)', dvhparam_label).group(1)            # string between brackets ()
        untangle_status = 1
    elif dvhparam_label[0:1] in ['V', 'D', 'C','H']:
        output_quantity = re.findall('([a-zA-Z]*)\d*.*', dvhparam_label)[0]     # string till numeric
        input_value = re.findall(r"[-+]?\d*\.?\d+", dvhparam_label)[0]              # old: r'\d+'     does not work for decimals
        if output_quantity == 'HI':
            input_unit = re.search(input_value + '(.*)', dvhparam_label).group(1)  # string after input_value (no space at end)
        else:
            input_unit = re.search(input_value + '(.*) ', dvhparam_label).group(1)      # string after input_value and before space
            output_unit = re.search('\((.*)\)', dvhparam_label).group(1)            # string between brackets ()
        input_value = float(input_value)
        untangle_status = 1

    # print
    if verbose:
        if untangle_status:
            print(f"\t\t\toutput_quantity: {output_quantity}   input_value: {input_value}   input_unit: {input_unit}   output_unit: {output_unit}")
        else:
            print(f"\t\t\tUnknown format for dvhparam_label")


    return(output_quantity, output_unit, input_value, input_unit, untangle_status)

def evaluatePlans(pathDVHParams, planNames, planTimes={}, planValues={}, planExtras={}, pathPatient=None, contourType=config.KEYNAME_CONTOUR_CLINICAL, save=True, verbose=False):

    # Step 0 - Initialize
    res = {}

    try:
        if Path(pathDVHParams).exists():
            
            # Step 1 - Get metrics (and other common vars)
            rois_and_dvhparams = read_dvhparamlist_csv_or_txt(str(pathDVHParams), ';')
            patient    = rayStationSave()
            case       = patient.Cases[0]
            exam       = case.Examinations[0]
            rois_case = [roi.Name for roi in case.PatientModel.RegionsOfInterest]
            res = {
                config.KEYNAME_PLANS: {planName : {} for planName in planNames}
                , config.KEYNAME_PLAN_DEBUGINFO: {planName: {
                    config.KEYNAME_TIME: planTimes.get(planName, -1)
                    , config.KEYNAME_OBJ_VALUE: planValues.get(planName, -1)
                } for planName in planNames}
                , config.KEYNAME_PLAN_EXTRAS: planExtras
            }
            plans = [ case.TreatmentPlans[planName] if checkForRTPlan(case, planName) else None for planName in planNames  ] # check if plan exists 

            # Step 2 - Loop over plans
            pt_all_plan_all_dvhparams_values = []
            print ('\n\n =================================================================\n\n')
            print (' - pathPatient: ', pathPatient)
            for planId, plan in enumerate(plans):
                try:
                    
                    # Step 2.1 - Check if plan exists
                    if plan is None:
                        print (' - [evaluatePlans()] Plan: {} does not exist'.format(planNames[planId]))
                        continue
                    print ('\n - [evaluatePlans()] --------------------- Plan: {} \n'.format(plan.Name))

                    # Step 2.2 - Check if dose exists
                    beamset = plan.BeamSets[plan.Name]
                    validDose = checkForRTDose(beamset)
                    if validDose:
                        
                        dose = plan.TreatmentCourse.TotalDose
                        dose.UpdateDoseGridStructures()
                        doselevels_processed = [beamset.Prescription.DosePrescriptions[0].DoseValue]
                        doselevels_processed =  [int(value) for dl, value in re.findall("(DL[0-9]) ([0-9]+)", beamset.Comment)]

                        # Step 3 - Loop over metrics
                        rois = []                               # holds per patient/plan all roi names
                        rois_synonym_hits = []                  # holds per patient/plan the roi synonym hit (if roi synonym input was used, and not the first roi name entry was found but a later entry)
                        headerline = []                         # holds the output headerline (ie all roi+dvh parameter labels) (headerline is filled every pt/plan iteration, only written once when output file is first generated)
                        pt_curplan_all_dvhparams_labels = []
                        pt_curplan_all_dvhparams_values = []    # holds per patient all dvh parameters of all the rois (an entire row in the output file)
                        for i, (roi_input, dvhparams_labels) in enumerate(rois_and_dvhparams.items()):
                            res[config.KEYNAME_PLANS][plan.Name][roi_input] = {}
                            if verbose: 
                                print(f"\n\tentry ({i+1}/{len(rois_and_dvhparams)}): {roi_input}")

                            try:
                            # first check for prescribeddose, or doselevels, if not then assume it is a roi (or roi synonym list)
                                if roi_input == 'prescribeddose#':
                                    prescribed_dose_value = beamset.Prescription.DosePrescriptions[0].DoseValue # assume one beamset, assume one doseprescription
                                    rois.append('Prescribed Dose')
                                    rois_synonym_hits.append('')
                                    pt_curplan_all_dvhparams_labels.append('')
                                    headerline.append('Prescribed Dose')
                                    pt_curplan_all_dvhparams_values.append(str(prescribed_dose_value))
                                elif roi_input == 'doselevels#':
                                    rois.append('Doselevels')
                                    rois_synonym_hits.append('')
                                    pt_curplan_all_dvhparams_labels.append('')
                                    headerline.append('Doselevels')
                                    pt_curplan_all_dvhparams_values.append(', '.join([str(x) for x in doselevels_processed]))
                                else:
                                    # we are dealing with a ROI (by exclusion)
                                    # process synonyms
                                    roi_input = roi_input.replace(', ',',')     # 'PTV_DL1_DVH,PTV_DL1'
                                    roi_synonym_list = roi_input.split(',')     # ['PTV_DL1_DVH', 'PTV_DL1']
                                    roi_synonym_master = roi_synonym_list[0]    # 'PTV_DL1_DVH'
                                    roi_synonym_hit = ''

                                    # check if current roi exists
                                    roi_status = 'inexistent'
                                    for roi_synonym in roi_synonym_list:
                                        if roi_synonym not in rois_case:
                                            # roi synonym not found, continue in loop
                                            if verbose: print(f'\t\t{roi_synonym}: {roi_status}')
                                        elif not case.PatientModel.StructureSets[exam.Name].RoiGeometries[roi_synonym].HasContours():
                                            # roi synonym hit WITHOUT contours, continue in loop
                                            roi_status = 'empty'
                                            roi_synonym_hit = roi_synonym
                                            if verbose: print(f'\t\t{roi_synonym}: {roi_status}')
                                        else:
                                            # roi synonym hit AND contoured, step out
                                            roi_status = 'contoured'
                                            roi_synonym_hit = roi_synonym
                                            if verbose: print(f'\t\t{roi_synonym}: {roi_status}')
                                            break
                                    
                                    # loop over the dvh parameters of this roi
                                    for dvhparam_label in dvhparams_labels:
                                        res[config.KEYNAME_PLANS][plan.Name][roi_input][dvhparam_label] = {}
                                        if not validDose:
                                            dvhparam_value = 'no dose'
                                        elif roi_status == 'contoured':
                                            dvhparam_value  = process_dvhparam(plan, doselevels_processed, dose, roi_synonym_hit, dvhparam_label, verbose=verbose)
                                        else:
                                            dvhparam_value = roi_status
                                        res[config.KEYNAME_PLANS][plan.Name][roi_input][dvhparam_label] = dvhparam_value

                                        # append roi status to list
                                        rois.append(roi_synonym_master)
                                        if roi_synonym is not roi_synonym_master:
                                            # only store synonyms (not the roi_synonym_master (first entry in roi_synonym_list))
                                            rois_synonym_hits.append(roi_synonym_hit)
                                        else:
                                            rois_synonym_hits.append('')
                                        pt_curplan_all_dvhparams_labels.append(dvhparam_label)
                                        headerline.append(roi_synonym_master + ' ' + dvhparam_label)
                                        pt_curplan_all_dvhparams_values.append(str(dvhparam_value))
                            
                            except:
                                traceback.print_exc()
                                print (' - [evaluatePlans()] Error in roi: ', roi_input)
                                

                        # Step 3.99 - Save plan information
                        pt_all_plan_all_dvhparams_values.append(pt_curplan_all_dvhparams_values)

                    else:
                        print (' - [evaluatePlans()] No RTDose found for plan: ', plan.Name)

                except:
                    traceback.print_exc()
                    print (' - [evaluatePlans()] Error in plan: ', plan.Name)
                    

            # Step 5 - Save
            if save:
                try:
                    import json
                    if pathPatient is not None:
                        if contourType == config.KEYNAME_CONTOUR_CLINICAL:
                            pathPlanStats = Path(pathPatient, config.KEY_PLAN_STATS_JSON)
                        elif contourType == config.KEYNAME_CONTOUR_AUTO:
                            pathPlanStats = Path(pathPatient, config.KEY_PLAN_STATS_JSON_AUTO)
                        elif contourType == config.KEYNAME_CONTOUR_ALL:
                            pathPlanStats = Path(pathPatient, config.KEY_PLAN_STATS_JSON_ALL)
                            
                        with open(str(pathPlanStats), 'w') as fp:
                            json.dump(res, fp, indent=4)
                except:
                    traceback.print_exc()

            # Step 6 - Print
            res2 = {}
            for planName in res[config.KEYNAME_PLANS]:
                res2[planName] = {}
                for roiName in res[config.KEYNAME_PLANS][planName]:
                    for dvhParam in res[config.KEYNAME_PLANS][planName][roiName]:
                        res2[planName][roiName + ' ({})'.format(dvhParam)] = res[config.KEYNAME_PLANS][planName][roiName][dvhParam]
            
            import pandas as pd
            df = pd.DataFrame.from_dict(res2, orient='index').T
            print (df)

        else:
            print (' - [evaluatePlans()] DVH params file does not exist: ', pathDVHParams)

    except:
        traceback.print_exc()
        

    return res

# rs_objective_template/helpers/condition_validator.py
class ConditionalElementValidator:
    
    def __init__(self, beamset, manual_selection_mapping=None):
        """ Helper class to check elements for conditions.
        Returns the element if no condition or condition is met, None otherwise

        Args:
            beamset: RayStation beamset object
        """
        self.beamset = beamset
        self.manual_selection_mapping = manual_selection_mapping

    def is_valid(self, element_tree: ElementTree):
        """ Checks ElementTree for conditions set in attributes

        Args:
            element_tree (Element): ElementTree to be checked

        Returns:
            bool: True if ElementTree passes conditions
        """
        validations = [True]

        conditional_fraction_doses = element_tree.get("conditionalFractionDoses")
        if conditional_fraction_doses:
            validations.append(self._evaluate_fraction_dose(conditional_fraction_doses))

        conditional_volume_roi = element_tree.get("conditionalVolumeRoi")
        if conditional_volume_roi:
            validations.append(self._evaluate_roi_volume(conditional_volume_roi, element_tree.get("conditionalVolumeBetween")))

        conditional_doselevels = element_tree.get("conditionalDoseLevels")
        if conditional_doselevels:
            validations.append(self._evaluate_doselevels(conditional_doselevels))

        conditional_doselevel_count = element_tree.get("conditionalDoseLevelCount")
        if conditional_doselevel_count:
            validations.append(self._evaluate_doselevel_count(conditional_doselevel_count))

        return all(validations)

    def set_manual_selection_mapping(self, mapping: dict):
        self.manual_selection_mapping = mapping

    def _evaluate_fraction_dose(self, conditional_fraction_doses: str):
        """ Evaluates 'conditionalFractionDoses'.
        conditionalFractionDoses is one or a list of fraction doses separated by '|'.

        Args:
            conditional_fraction_doses (str): Conditional fraction doses string separated by '|'

        Returns:
            bool: True if check passes
        """
        fraction_dose = round(self.beamset.Prescription.PrimaryDosePrescription.DoseValue / self.beamset.FractionationPattern.NumberOfFractions, ndigits=1)
        valid_fraction_doses = [float(fx) for fx in conditional_fraction_doses.split("|")]
        valid = fraction_dose in valid_fraction_doses
        print(f"   ** [ConditionalElementValidator._evaluate_fraction_dose()] Fraction dose condition for {fraction_dose} {'meets' if valid else 'does not meet'} condition {valid_fraction_doses}")
        return valid

    def _evaluate_roi_volume(self, conditional_roi: str, conditional_volume_range: str):
        """ Evaluates 'conditionalVolumeRoi'.
        conditionalVolumeRoi is a string with the name of the ROI 'conditionalVolumeBetween' compares the volume to.

        Args:
            conditional_roi (str):  Roi to compare volume for
            conditional_volume_range (str):  Volume range string to compare, seperated by '-'

        Returns:
            bool: True if check passes
        """
        if self.manual_selection_mapping:
            conditional_roi = self.manual_selection_mapping.get(conditional_roi, None)

        conditional_volume = [float(v) for v in conditional_volume_range.split("-")]

        structure_set = self.beamset.GetStructureSet()
        if conditional_roi not in [roi.OfRoi.Name for roi in structure_set.RoiGeometries]:
            return False
        roi_volume = round(structure_set.RoiGeometries[conditional_roi].GetRoiVolume(), ndigits=2)
        valid = conditional_volume[0] <= roi_volume < conditional_volume[1]
        print(f"Roi volume {conditional_roi}({roi_volume}cc) {'meets' if valid else 'does not meet'} condition 'between {conditional_volume}cc'")
        return valid

    def _evaluate_doselevels(self, conditional_doselevels: str):
        """ Evaluates 'conditionalDoseLevels'.
        conditionalDoseLevels is one or a list of doselevels separated by '|'.

        Args:
            conditional_doselevels (str): Conditional doselevels string separated by '|'

        Returns:
            bool: True if check passes
        """
        doselevels = {dl: float(value) for dl, value in re.findall("(DL[0-9]) ([0-9]+)", self.beamset.Comment)}
        valid_doselevels = conditional_doselevels.split("|")
        valid = any([f"{k} {int(v)}" in valid_doselevels for k, v in doselevels.items()])
        print(f"   ** [ConditionalElementValidator._evaluate_doselevels()] Doselevel condition for {doselevels} {'meets' if valid else 'does not meet'} condition {valid_doselevels}")
        return valid

    def _evaluate_doselevel_count(self, doselevel_count: int):
        """ Evaluates 'conditionalDoseLevelCount'.
        conditionalDoseLevelCount is an integer defining how many doselevels must be present.

        Args:
            doselevel_count (int): Conditional doselevel count

        Returns:
            bool: True if check passes
        """
        doselevels = {dl: float(value) for dl, value in re.findall("(DL[0-9]) ([0-9]+)", self.beamset.Comment)}
        valid = len(doselevels.keys()) == int(doselevel_count)
        print(f"   -- [INFO][ConditionalElementValidator._evaluate_doselevel_count()] Doselevel count condition for {doselevels} {'meets' if valid else 'does not meet'} condition {doselevel_count}")
        return valid

# rs_objective_template/models/objective.py
class InvalidObjectiveException(Exception):
    pass

class InvalidDoseLevelException(Exception):
    pass

class Objective(object):
    def __new__(cls, plan, beamset, roi: str, optimization_tree: ElementTree, condition_validator: ConditionalElementValidator = None):
        """ Objective Factory, returns different Objective object depending on functionType of optimization_tree

        Args:
            plan: RayStation plan
            beamset: RayStation beamset
            roi (str): ROI name to apply objective to
            optimization_tree (Element): Optimization XML tree
        """
        function_type = optimization_tree.get("functionType")
        if function_type == "DoseFallOff":
            return ObjectiveDoseFallOff(plan, beamset, roi, optimization_tree, condition_validator)
        if function_type in ["MinDose", "MaxDose", "UniformDose"]:
            return ObjectiveDose(plan, beamset, roi, optimization_tree, condition_validator)
        if function_type in ["MinDvh", "MaxDvh"]:
            return ObjectiveDoseVolume(plan, beamset, roi, optimization_tree, condition_validator)
        if function_type in ["MaxEud", "MinEud", "UniformEud"]:
            return ObjectiveEud(plan, beamset, roi, optimization_tree, condition_validator)
        if function_type == "UniformityConstraint":
            return ObjectiveUniformity(plan, beamset, roi, optimization_tree, condition_validator)

        raise InvalidObjectiveException(f"'{function_type}' is not a valid FunctionType")

class ObjectiveBase:
    def __init__(self, plan, beamset, roi: str, optimization_tree: ElementTree, condition_validator: ConditionalElementValidator = None):
        """ Objective base class, provides all base functionality for Objectives

        Args:
            plan: RayStation plan
            beamset: RayStation beamset
            roi (str): ROI name to apply objective to
            optimization_tree (Element): Optimization XML tree
        """
        self.plan = plan
        self.beamset = beamset
        self.ConditionValidator = condition_validator or ConditionalElementValidator(plan, beamset)
        # Parse doselevels in beamset comment and put in dict for easy lookup
        self.doselevels = {dl: float(value) for dl, value in re.findall("(DL[0-9]) ([0-9]+)", self.beamset.Comment)}
        self.roi_name = roi
        self.optimization_tree = optimization_tree
        self.function_type = self.optimization_tree.get("functionType")
        self.is_constraint = self._get_bool(self.optimization_tree.get("isConstraint"))
        self.is_robust = self._get_bool(self.optimization_tree.get("isRobust"))

        if self.optimization_tree.get("restrictToBeams") and self.beamset.Modality != "Protons":
            raise InvalidObjectiveException("'restrictToBeams' only valid for proton plans...")
        self.restrict_to_beams = [
            beam.Name for beam in self.beamset.Beams \
                if beam.Name in self.optimization_tree.get("restrictToBeams").split(',')
                ] if self.optimization_tree.get("restrictToBeams") not in [None, "All"] else self.optimization_tree.get("restrictToBeams") or False

        # if self.roi_name == 'Opt_CTV_L':
        #     print ('\n - [ObjectiveBase()] self.roi_name: ', self.roi_name)
        #     print (' - [ObjectiveBase()] self.restrict_to_beams: ', self.restrict_to_beams)
        #     pdb.set_trace()

        # Check parameter tree for conditions and get only valid tree or default, raises exception when None or multiple parameter trees are found
        self.parameters_tree = [fp for fp in self.optimization_tree.findall("FunctionParameters") if self.ConditionValidator.is_valid(fp)]
        if len(self.parameters_tree) == 2:
            self.parameters_tree = [pt for pt in self.parameters_tree if not pt.get("default")]
        if len(self.parameters_tree) == 1:
            self.parameters_tree = self.parameters_tree[0]
        else:
            raise InvalidObjectiveException(f"Multiple or no FunctionParameters found for OptimizationFunction '{self.function_type}' in '{self.roi_name}'...")

        if not self.parameters_tree.get("weight"):
            raise InvalidObjectiveException(f"No weight specified for '{self.function_type}' in '{self.roi_name}'...")
        self.weight = float(self.parameters_tree.get("weight"))

    def _get_doselevel(self, doselevel_tree_name: str):
        """ Parses a doselevel tree by name of element (i.e. DoseLevel, LowDoseLevel, HighDoseLevel)
        When relativeToDoseLevel attribute is specified, calculates the absolute dose of the objective for the specified dose level or 'Highest'

        Args:
            doselevel_tree_name (str): Name of element to be parsed

        Returns:
            int: Absolute dose of doselevel element in cGy
        """
        doselevel_tree = [dl for dl in self.parameters_tree.findall(doselevel_tree_name) if self.ConditionValidator.is_valid(dl)]
        if len(doselevel_tree) == 2:
            doselevel_tree = [dl for dl in doselevel_tree if not dl.get("default")]
        if len(doselevel_tree) == 1:
            doselevel_tree = doselevel_tree[0]
        else:
            raise InvalidObjectiveException(f"None or multiple DoseLevels found for FunctionParameter '{self.function_type}' in '{self.roi_name}'...")

        relative_to_doselevel = doselevel_tree.get("relativeToDoseLevel", default=False)
        if not relative_to_doselevel:
            # return int(doselevel_tree.text)
            return float(doselevel_tree.text)

        if relative_to_doselevel.lower() == "highest":
            relative_to_doselevel = max(self.doselevels, key=self.doselevels.get)
        try:
            relative_total_dose = self.doselevels[relative_to_doselevel]
        except KeyError:
            raise InvalidDoseLevelException(f"Relative doselevel '{relative_to_doselevel}' not found!")

        return round(relative_total_dose * (float(doselevel_tree.text)/100))

    def _get_bool(self, boolstring: str):
        """ Parses string to boolean value (wrapper around strtobool), returns False when None

        Args:
            boolstring (str): "true", "false", "yes", "no", etc

        Returns:
            bool
        """
        return bool(strtobool(boolstring)) if boolstring else False

    def _create_optimization_function(self, restrict_to_beam=None):
        """ Creates the RayStation optimization function and returns the reference to that function

        Returns:
            OptimizationFunction: RayStation optimization function object
        """
        # Set default values
        restrict_to_beamset = None
        restrict_all_beams_individually = False

        if self.restrict_to_beams == "All":
            restrict_to_beamset = self.beamset.DicomPlanLabel
            restrict_all_beams_individually = True

        if restrict_to_beam:
            n = self.beamset.DicomPlanLabel

        # PlanOptimizations index is related to BeamSet index
        return self.plan.PlanOptimizations[self.plan.BeamSets.IndexOf(self.beamset)].AddOptimizationFunction(
            FunctionType=self.function_type
            , RoiName=self.roi_name
            , IsConstraint=self.is_constraint
            , IsRobust=self.is_robust
            , RestrictToBeam=restrict_to_beam
            , RestrictToBeamSet=restrict_to_beamset
            , RestrictAllBeamsIndividually=restrict_all_beams_individually)

    def _set_function_parameters(self, optimization_function):
        """ Sets base function parameters of specified optimization function

        Args:
            optimization_function: RayStation optimization function
        """
        optimization_function.DoseFunctionParameters.Weight = self.weight
        # print ('weight', self.weight)

    # def apply(self):
    #     """ Start _create_optimization_function and passes reference to _set_function_parameters """
    #     optimization_function = self._create_optimization_function()
    #     self._set_function_parameters(optimization_function)
    
    def apply(self):
        """ Start _create_optimization_function and passes reference to _set_function_parameters """
        self.optimization_function = self._create_optimization_function()
        self._set_function_parameters(self.optimization_function)
    
    def update_weight(self, weight):
        self.weight = weight
        self.optimization_function.DoseFunctionParameters.Weight = weight

class ObjectiveDoseFallOff(ObjectiveBase):
    def __init__(self, *args, **kwargs):
        """ DoseFallOff objective, adds DoseFallOff specific attributes"""
        super().__init__(*args, **kwargs)
        self.adapt_to_target_doselevels = self._get_bool(self.parameters_tree.get("adaptToTargetDoseLevels"))
        self.low_dose_distance = float(self.parameters_tree.get("lowDoseDistance"))
        self.high_doselevel = self._get_doselevel("HighDoseLevel")
        self.low_doselevel = self._get_doselevel("LowDoseLevel")

    def _set_function_parameters(self, optimization_function):
        """ Override to add DoseFallOff specific attributes

        Args:
            optimization_function: RayStation optimization function
        """
        super()._set_function_parameters(optimization_function)
        optimization_function.DoseFunctionParameters.AdaptToTargetDoseLevels = self.adapt_to_target_doselevels
        optimization_function.DoseFunctionParameters.LowDoseDistance = self.low_dose_distance
        optimization_function.DoseFunctionParameters.HighDoseLevel = self.high_doselevel # in cGy
        optimization_function.DoseFunctionParameters.LowDoseLevel = self.low_doselevel # in cGy

class ObjectiveDose(ObjectiveBase):
    def __init__(self, *args, **kwargs):
        """ Dose objective, adds Dose specific attributes"""
        super().__init__(*args, **kwargs)
        self.doselevel = self._get_doselevel("DoseLevel")

    def _set_function_parameters(self, optimization_function):
        """ Override to add Dose specific attributes

        Args:
            optimization_function: RayStation optimization function
        """
        super()._set_function_parameters(optimization_function)
        optimization_function.DoseFunctionParameters.DoseLevel = self.doselevel

    def apply(self):
        """ Override for applying multiple objectives when `restrictToBeams` is set
        Divides the doselevel of the objective by the number of beams in `restrictToBeams`
        """
        
        # Just run the default apply action if not restricted to individual beams
        if not self.restrict_to_beams or self.restrict_to_beams == "All":
            super().apply()
            return

        # Divide doselevel of function by number of beams function is restricted to. Create a objective function per restricted beam
        self.doselevel = self.doselevel / len(self.restrict_to_beams)
        for beam in self.restrict_to_beams:
            
            optimization_function = self._create_optimization_function(restrict_to_beam=beam)
            self._set_function_parameters(optimization_function)

class ObjectiveDoseVolume(ObjectiveDose):
    def __init__(self, *args, **kwargs):
        """ Dose objective, adds DoseVolume specific attributes"""
        super().__init__(*args, **kwargs)
        self.percent_volume = float(self.parameters_tree.get("percentVolume"))

    def _set_function_parameters(self, optimization_function):
        """ Override to add DoseVolume specific attributes

        Args:
            optimization_function: RayStation optimization function
        """
        super()._set_function_parameters(optimization_function)
        optimization_function.DoseFunctionParameters.PercentVolume = self.percent_volume

class ObjectiveEud(ObjectiveDose):
    def __init__(self, *args, **kwargs):
        """ Dose objective, adds Eud specific attributes"""
        super().__init__(*args, **kwargs)
        self.eud_parameter_a = float(self.parameters_tree.get("eudParameterA"))

    def _set_function_parameters(self, optimization_function):
        """ Override to add Eud specific attributes

        Args:
            optimization_function: RayStation optimization function
        """
        super()._set_function_parameters(optimization_function)
        optimization_function.DoseFunctionParameters.EudParameterA = self.eud_parameter_a

class ObjectiveUniformity(ObjectiveBase):
    def __init__(self, *args, **kwargs):
        """ Dose objective, adds Uniformity specific attributes"""
        super().__init__(*args, **kwargs)
        self.percent_std_deviation = float(self.parameters_tree.get("percentStdDeviation"))

    def _set_function_parameters(self, optimization_function):
        """ Override to add Uniformity specific attributes

        Args:
            optimization_function: RayStation optimization function
        """
        super()._set_function_parameters(optimization_function)
        optimization_function.DoseFunctionParameters.PercentStdDeviation = self.percent_std_deviation

# [entry point] rs_objective_template/helpers/objective_template_manager.py
class ObjectiveTemplateManager:

    def __init__(self, plan, beamset):
        """ Manager object for Objective templates, needs RayStation plan and beamset objects to get ROIs
         and plan information for validation of conditions

        Args:
            plan: RayStation plan
            beamset: RayStation beamset
        """

        self.plan    = plan
        self.beamset = beamset
        self.all_rois_in_case = [roi.OfRoi.Name for roi in self.beamset.GetStructureSet().RoiGeometries]
        self.rois_in_case = [roi.OfRoi.Name for roi in self.beamset.GetStructureSet().RoiGeometries if roi.HasContours()]
        
        # Instantiate ConditionValidator with beamset for validation of conditions (i.e fraction dose, roi volume)
        self.ConditionValidator = ConditionalElementValidator(self.beamset)
        self.objectives = []
        self.manual_roi_mapping = {}
        self.template_name = None

    def parse_xml(self, path_to_xml):
        """ Parses XML file into list of Objective objects

        Args:
            path_to_xml (str): Absolute path to XML template file
        """
        self.objectives.clear()
        self.template_name = os.path.basename(path_to_xml)[:-4]
        objective_tree = ElementTree.parse(path_to_xml).find("ObjectiveTemplate")
        valid_rois_tree = []
        manual_rois_tree = []
        unmatched_rois_tree = []

        # Loop over ROI's and validate conditions
        for roi in [roi for roi in objective_tree.findall("Roi") if self.ConditionValidator.is_valid(roi)]:
            roi_name = roi.get("name")
            manual_roi = roi.get("userSelection")
            if manual_roi:
                manual_rois_tree.append(roi)
            elif roi_name not in self.all_rois_in_case:
                unmatched_rois_tree.append(roi)
            elif not roi_name:
                continue
            else:
                valid_rois_tree.append(roi)

        # If manual ROI's or unmatched ROI's; display selection dialog window and add selection to valid ROI list
        if 0:
            if manual_rois_tree or unmatched_rois_tree:
                with RoiSelectDialog(self.all_rois_in_case) as roi_select_dialog:
                    if roi_select_dialog.ShowDialog(manual_rois_tree, unmatched_rois_tree, self.template_name):
                        valid_rois_tree[0:0] = [tree for tree in roi_select_dialog.assigned_roi_tree_list if tree.get("name") != "--"]
                        self.manual_roi_mapping = {tree.get('userSelection'): tree.get("name") for tree in roi_select_dialog.assigned_roi_tree_list}
                        self.ConditionValidator.set_manual_selection_mapping(self.manual_roi_mapping)
                    else:
                        return
        else:
            self.manual_roi_mapping = {}
            for roi in objective_tree.findall("Roi"):
                self.manual_roi_mapping[roi] = roi 
            # print (self.manual_roi_mapping)
            self.ConditionValidator.set_manual_selection_mapping(self.manual_roi_mapping)

        # Add all valid objectives in valid ROI's to objectives list
        for roi in valid_rois_tree:
            for objective_tree in [ot for ot in roi.findall("OptimizationFunction") if self.ConditionValidator.is_valid(ot)]:
                objective = Objective(self.plan, self.beamset, roi.get("name"), objective_tree, self.ConditionValidator)
                self.objectives.append(objective)

# rs_isodose_template
class ConditionalElementValidatorForIsoDose:
    def __init__(self, doselevels: dict):
        """ Helper class to check elements for conditions.
        Returns the element if no condition or condition is met, None otherwise

        Args:
            doselevels: Dict of doselevels in format {"DL1": 3500, "DL2": 4000}
        """
        self.doselevels = doselevels

    def is_valid(self, element_tree: ElementTree):
        """ Checks ElementTree for conditions set in attributes

        Args:
            element_tree (Element): ElementTree to be checked

        Returns:
            bool: True if ElementTree passes conditions
        """
        if element_tree.get("conditionalDoseLevels"):
            return self._evaluate_doselevels(element_tree)

        return True

    def _evaluate_doselevels(self, element_tree: ElementTree):
        """ Evaluates 'conditionalDoseLevels'.
        conditionalDoseLevels is a string with the name of the DoseLevels present (separated by ,) i.e. "DL1, DL2, DL3".

        Args:
            element_tree (Element):  ElementTree to be checked

        Returns:
            bool: True if check passes
        """
        conditional_doselevels = element_tree.get("conditionalDoseLevels")
        if conditional_doselevels:
            if len(self.doselevels) != len(conditional_doselevels.split(",")):
                return False
            return all([dl in conditional_doselevels.split(",") for dl in self.doselevels.keys()])

        return True

class IsodoseManager:
    def __init__(self, path_isodose_xml, case, doselevels: dict, template_xml_path=None):
        self.case = case
        self.doselevels = doselevels
        self.max_doselevel = self.doselevels[max(self.doselevels, key=lambda k: self.doselevels[k])]
        self.ConditionValidator = ConditionalElementValidatorForIsoDose(doselevels)
        self.isodose_tree = []

        if Path(path_isodose_xml).exists():
            self._parse_xml(str(path_isodose_xml))

    def _parse_xml(self, path_to_xml: str):
        """ Parses XML file into dict of ColorMaps
        Args:
            path_to_xml (str): Absolute path to XML template file
        """
        isodose_trees = [template for template in ElementTree.parse(path_to_xml).findall("IsodoseTemplate") if self.ConditionValidator.is_valid(template)]

        for tree in isodose_trees:
            self.isodose_tree.extend(tree.findall("Isodose"))

    def map_isodose(self):
        isodose_dict = {}
        for isodose_element in self.isodose_tree:
            isodose = Isodose(isodose_element, self.doselevels)
            isodose_dict.update(isodose.isodose_dict)

        self.case.CaseSettings.DoseColorMap.PresentationType = "Absolute"
        self.case.CaseSettings.DoseColorMap.ReferenceValue = self.max_doselevel
        try:
            self.case.CaseSettings.DoseColorMap.ColorMapReferenceType = "RelativePrescription"
        except:
            self.case.CaseSettings.DoseColorMap.ColorMapReferenceType = "ReferenceValue"
        self.case.CaseSettings.DoseColorMap.ColorTable = isodose_dict

class Isodose:
    def __init__(self, isodose_element: ElementTree, doselevels: dict):
        self.color = [int(color) for color in isodose_element.get("color").split(",")]
        percentage = int(isodose_element.text)
        relative_doselevel = isodose_element.get("relativeToDoseLevel")
        absolute_dose = round(doselevels[relative_doselevel] * (percentage / 100))
        self.relative_percentage = (absolute_dose / doselevels[max(doselevels, key=lambda k: doselevels[k])]) * 100

    @property
    def isodose_dict(self):
        from System.Drawing import Color # is importable within raystation
        return {self.relative_percentage: Color.FromArgb(*self.color)}

##########################################################################################
#                                     MICHAELS CODE                                      #
##########################################################################################

# Import python class libs
from enum import Enum
from typing import Dict
from abc import ABC, abstractmethod

# rs_ntcp_kno/models/enums.py
class TreatmentType(Enum):
    PRIMARY = "Primary"
    POSTOPERATIVE = "Postoperative"

# rs_ntcp_kno/models/model_roi.py
class ModelRoi:
    REMOVED_NAME = "--Removed--"

    def __init__(self, name: str = None, color: Dict[str, int] = None):
        self._name = name
        self.dose = None
        self.volume = None
        self.color = color

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if self._name == value:
            return
        self._name = value
        self.dose = None
        self.volume = None

# rs_ntcp_kno/models/ntcp_base.py
class InvalidRoiDoseException(Exception):
    ...

# rs_ntcp_kno/models/ntcp_base.py
class NtcpKnoAbstractBase(ABC):
    """ NTCP model interface (abstract base class)"""

    class ModelRoiType(Enum):
        """ RoiTypes needed for the model
        Model specific and must be overridden in subclass
        """
        ...

    MODEL_ROI_COLOR_MAPPING: Dict[ModelRoiType, Dict[str, int]] = None

    class BaselineType(Enum):
        """ Baseline types for the model
        Model specific and must be overridden in subclass
        """
        ...

    @abstractmethod
    def __init__(self, plan, grade: int):
        self.model_name = None
        self.ntcp_proton_criterion = None
        self.plan = plan
        self._grade = grade
        self._baseline_type = None
        self._treatment_type = None
        for roi_type in self.ModelRoiType:
            setattr(self, roi_type.name.lower(), ModelRoi(color=self.MODEL_ROI_COLOR_MAPPING[roi_type]))

    @property
    def model_rois(self) -> Dict[ModelRoiType, ModelRoi]:
        """ Property containing a dict of model ROI's with a bound ModelRoi instance
        Returns:
            ModelRoiType mappings to ModelRoi instances ie. {ModelRoiType: ModelRoi}
        """
        roi_dict = {}
        for roi_type in self.ModelRoiType:
            roi_dict[roi_type] = getattr(self, roi_type.name.lower())
        return roi_dict

    def set_model_roi(self, roi_name: str, roi_type: ModelRoiType):
        """ Bind a ROI name to a model ROI type
        Args:
            roi_name: Name of the ROI to bind
            roi_type: ModelRoiType of this model to bind to
        """
        if not roi_name:
            return
        for key, value in self.model_rois.items():
            if key != roi_type and value.name == roi_name:
                self.model_rois[roi_type].name = None
                return
        self.model_rois[roi_type].name = roi_name

    def set_baseline(self, baseline: BaselineType):
        """ Set the baseline type of this model instance
        Args:
             baseline: BaselineType to use for this model instance
        """
        self._baseline_type = baseline

    def set_treatment_type(self, treatment_type: TreatmentType):
        """ Set the treatment type of this model instance
        Args:
            treatment_type: TreatmentType to use for this model instance
        """
        self._treatment_type = treatment_type

    @abstractmethod
    def get_ntcp(self, dose_overrides: Dict[ModelRoiType, float] = None) -> float:
        """ Get NTCP value of this model instance
        Args:
            dose_overrides: Dict of ModelRoiType mean dose overrides to use
        Returns:
            float: NTCP value
        """
        ...

    def get_ntcp_curve(self, roi_type: ModelRoiType, center_dose_gy: float = None, curve_size_gy: int = 20, step_size_gy: float = 0.1) -> Dict[float, float]:
        """ Get NCTP curve for specific ROI
        Args:
            roi_type: ModelRoiType to get NTCP curve for
            center_dose_gy: Dose in Gy to center curve around, defaults to mean dose of chosen ROI
            curve_size_gy: Range of the created curve in Gy
            step_size_gy: Step size in Gy of calculated points in the curve
        Returns:
            Dict of floats ie. {dose_in_gy: ntcp_value}
        """
        rounded_center_dose = round(center_dose_gy or self.get_mean_roi_dose(roi_type))
        dose_points = np.arange(max(rounded_center_dose - round(curve_size_gy / 2), 0), rounded_center_dose + round(curve_size_gy / 2) + 1, step_size_gy)
        ntcp_curve_values = {}

        for dose in dose_points:
            ntcp_curve_values[dose] = self.get_ntcp({roi_type: dose})

        return ntcp_curve_values

    def get_ntcp_per_gy(self, roi_type: ModelRoiType) -> float:
        """ Get NTCP gradient for ModelRoiType for current model instance
        Args:
            roi_type: ModelRoiType to get gradient for
        Returns:
            Gradient in NTCP/Gy for given ModelRoiType
        """
        point_before = self.get_ntcp({roi_type: max(self.get_mean_roi_dose(roi_type) - 0.5, 0)})
        point_after = self.get_ntcp({roi_type: max(self.get_mean_roi_dose(roi_type) + 0.5, 1)})
        return point_after - point_before

    def get_ntcp_threshold(self) -> float:
        """ Get NTCP PV threshold for current model instance
        Returns:
            NTCP threshold for plan comparison
        """
        return self.get_base_ntcp() + self.ntcp_proton_criterion

    def get_base_ntcp(self):
        """ Get NTCP at 0 Gy for current model instance
        Returns:
            NTCP value at 0 Gy
        """
        return self.get_ntcp({roi_type: 0 for roi_type, _ in self.model_rois.items()})

    # Helper functions
    def get_mean_roi_dose(self, roi_type: ModelRoiType) -> float:
        """ Get mean dose of ModelRoiType bound ROI
            Calculates only once, then stores it in bound model for future requests
        Args:
            roi_type: ModelRoiType to get mean dose for
        Returns:
            Mean dose of given ROI in Gy
        """
        roi_cache = self.model_rois[roi_type]
        if roi_cache.dose is None:
            if roi_cache.name == ModelRoi.REMOVED_NAME:
                roi_cache.dose = 0
                return roi_cache.dose
            roi_dose_grid = self.plan.TreatmentCourse.TotalDose.GetDoseGridRoi(RoiName=roi_cache.name)
            if roi_dose_grid is None:
                if roi_cache.name is None:
                    raise InvalidRoiDoseException(f"ROI could not be set, duplicate entry?")
                raise InvalidRoiDoseException(f"ROI '{roi_cache.name}' has no dose grid")
            if roi_dose_grid.RoiVolumeDistribution is None:
                raise InvalidRoiDoseException(f"ROI '{roi_cache.name}' has no volume distribution, update dose statistics?")
            roi_cache.dose = self.plan.TreatmentCourse.TotalDose.GetDoseStatistic(RoiName=roi_cache.name, DoseType="Average") / 100
        return roi_cache.dose

# rs_ntcp_kno/models/ntcp_xerostomia.py
class NtcpKnoXerostomia(NtcpKnoAbstractBase):
    PRIMARY_CONSTANTS: Dict[str, Dict[int, float]] = {
        "grade":         {2: -2.2951,   3: -3.7286},
        "parotid":       {2: 0.0996,    3: 0.0855},
        "submandibular": {2: 0.0182,    3: 0.0156},
        "xerostomia_0":  {2: 0,         3: 0},
        "xerostomia_1":  {2: 0.4950,    3: 0.4249},
        "xerostomia_2":  {2: 1.2070,    3: 1.0361}
    }
    POSTOPERATIVE_CONSTANTS: Dict[str, Dict[int, float]] = {
        "grade":         {2: -1.6824,   3: -4.3613},
        "parotid":       {2: 0.0388,    3: 0.1054},
        "submandibular": {2: 0.0071,    3: 0.0193},
        "xerostomia_0":  {2: 0,         3: 0},
        "xerostomia_1":  {2: 0.1925,    3: 0.5234},
        "xerostomia_2":  {2: 0.4695,    3: 1.2763}
    }

    class ModelRoiType(Enum):
        PAROTID_LEFT = "Parotid left"
        PAROTID_RIGHT = "Parotid right"
        SUBMANDIBULARS = "Submandibulars"

    MODEL_ROI_COLOR_MAPPING: Dict[ModelRoiType, Dict[str, int]] = {
        ModelRoiType.PAROTID_LEFT: {"r": 0, "g": 192, "b": 0},
        ModelRoiType.PAROTID_RIGHT: {"r": 0, "g": 192, "b": 192},
        ModelRoiType.SUBMANDIBULARS: {"r": 139, "g": 69, "b": 19}
    }

    class BaselineType(Enum):
        NONE = "None (helemaal niet)"
        LITTLE = "A little (een beetje)"
        SEVERE = "Severe (nogal/heel erg)"

    BASELINETYPE_MAPPING: Dict[BaselineType, int] = {
        BaselineType.NONE:   0,
        BaselineType.LITTLE: 1,
        BaselineType.SEVERE: 2
    }

    def __init__(self, plan, grade: int):
        super().__init__(plan, grade)
        self.model_name = f"NTCP Xerostomia Grade ≥ {grade}"
        self.ntcp_proton_criterion = 0.1 if grade == 2 else 0.05

    def get_ntcp(self, dose_overrides: Dict[ModelRoiType, float] = None) -> float:
        """ Get Xerostomia NTCP value

        Args:
            dose_overrides: Dict of ModelRoiType mean dose overrides to use
        Returns:
            float: NTCP value
        """
        if not all([self._treatment_type, self._baseline_type]):
            raise ValueError("Cannot calculate NTCP without setting all variables")

        constants = self.PRIMARY_CONSTANTS if self._treatment_type == TreatmentType.PRIMARY else self.POSTOPERATIVE_CONSTANTS

        intermediate_constant = constants["grade"][self._grade]
        intermediate_constant += constants[f"xerostomia_{self.BASELINETYPE_MAPPING[self._baseline_type]}"][self._grade]

        parotid_l_dose = dose_overrides.get(self.ModelRoiType.PAROTID_LEFT) if dose_overrides is not None and dose_overrides.get(self.ModelRoiType.PAROTID_LEFT) is not None else self.get_mean_roi_dose(self.ModelRoiType.PAROTID_LEFT)
        parotid_r_dose = dose_overrides.get(self.ModelRoiType.PAROTID_RIGHT) if dose_overrides is not None and dose_overrides.get(self.ModelRoiType.PAROTID_RIGHT) is not None else self.get_mean_roi_dose(self.ModelRoiType.PAROTID_RIGHT)
        submandibulars_dose = dose_overrides.get(self.ModelRoiType.SUBMANDIBULARS) if dose_overrides is not None and dose_overrides.get(self.ModelRoiType.SUBMANDIBULARS) is not None else self.get_mean_roi_dose(self.ModelRoiType.SUBMANDIBULARS)

        exponent = {"parotid_l":      constants["parotid"][self._grade] * math.sqrt(parotid_l_dose),
                    "parotid_r":      constants["parotid"][self._grade] * math.sqrt(parotid_r_dose),
                    "submandibulars": constants["submandibular"][self._grade] * submandibulars_dose}

        ntcp_total = 1 / (1 + math.exp(-intermediate_constant - sum(exponent.values())))

        return ntcp_total

# rs_ntcp_kno/models/ntcp_dysphagia.py
class NtcpKnoDysphagia(NtcpKnoAbstractBase):
    PRIMARY_CONSTANTS: Dict[str, Dict[int, float]] = {
        "grade":                    {2: -4.0536,    3: -7.6174},
        "oral_cavity":              {2: 0.03,       3: 0.0259},
        "pcm_superior":             {2: 0.0236,     3: 0.0203},
        "pcm_middle":               {2: 0.0095,     3: 0.0303},
        "pcm_inferior":             {2: 0.0133,     3: 0.0341},
        "dysphagia_0":              {2: 0,          3: 0},
        "dysphagia_1":              {2: 0.9382,     3: 0.5738},
        "dysphagia_2":              {2: 1.29,       3: 1.4718},
        "tumor_loc_oral_cavity":    {2: 0,          3: 0},
        "tumor_loc_pharynx":        {2: -0.6281,    3: 0.0387},
        "tumor_loc_larynx":         {2: -0.7711,    3: -0.5303}
    }

    POSTOPERATIVE_CONSTANTS: Dict[str, Dict[int, float]] = {
        "grade":                    {2: -2.4138,    3: -3.2594},
        "oral_cavity":              {2: 0.0192,     3: 0.0063},
        "pcm_superior":             {2: 0.0151,     3: 0.0050},
        "pcm_middle":               {2: 0.0060,     3: 0.0074},
        "pcm_inferior":             {2: 0.0085,     3: 0.0084},
        "dysphagia_0":              {2: 0,          3: 0},
        "dysphagia_1":              {2: 0.5985,     3: 0.1404},
        "dysphagia_2":              {2: 0.8227,     3: 0.3603},
        "tumor_loc_oral_cavity":    {2: 0,          3: 0},
        "tumor_loc_pharynx":        {2: -0.4007,    3: 0.0095},
        "tumor_loc_larynx":         {2: -0.4918,    3: -0.1298}
    }

    class ModelRoiType(Enum):
        ORAL_CAVITY = "Oral cavity"
        PCM_SUPERIOR = "Superior PCM"
        PCM_MIDDLE = "Middle PCM"
        PCM_INFERIOR = "Inferior PCM"

    MODEL_ROI_COLOR_MAPPING: Dict[ModelRoiType, Dict[str, int]] = {
        ModelRoiType.ORAL_CAVITY: {"r": 192, "g": 192, "b": 0},
        ModelRoiType.PCM_SUPERIOR:  {"r": 114, "g": 114, "b": 255},
        ModelRoiType.PCM_MIDDLE: {"r": 138, "g": 224, "b": 52},
        ModelRoiType.PCM_INFERIOR: {"r": 204, "g": 100, "b": 204}
    }

    class BaselineType(Enum):
        GRADE_0_1 = "Grade 0-1"
        GRADE_2 = "Grade 2"
        GRADE_3_4 = "Grade 3-4"

    BASELINETYPE_MAPPING = {
        BaselineType.GRADE_0_1:   0,
        BaselineType.GRADE_2: 1,
        BaselineType.GRADE_3_4: 2
    }

    class TumorLocationType(Enum):
        ORAL_CAVITY = "Oral cavity"
        PHARYNX = "Pharynx"
        LARYNX = "Larynx"

    def __init__(self, plan, grade: int):
        super().__init__(plan, grade)
        self.model_name = f"NTCP Dysphagia Grade ≥ {grade}"
        self.ntcp_proton_criterion = 0.1 if grade == 2 else 0.05
        self._tumor_location: NtcpKnoDysphagia.TumorLocationType = None

    def set_tumor_location(self, tumor_location: TumorLocationType):
        self._tumor_location = tumor_location

    def get_ntcp(self, dose_overrides: Dict[ModelRoiType, float] = None):
        """ Get Dysphagia NTCP value

        Args:
            dose_overrides: Dict of ModelRoiType mean dose overrides to use
        Returns:
            float: NTCP value
        """
        if not all([self._treatment_type, self._baseline_type, self._tumor_location]):
            raise ValueError("Cannot calculate NTCP without setting all variables")

        constants = self.PRIMARY_CONSTANTS if self._treatment_type == TreatmentType.PRIMARY else self.POSTOPERATIVE_CONSTANTS

        intermediate_constant = constants["grade"][self._grade]
        intermediate_constant += constants[f"dysphagia_{self.BASELINETYPE_MAPPING[self._baseline_type]}"][self._grade]
        intermediate_constant += constants[f"tumor_loc_{self._tumor_location.name.lower()}"][self._grade]

        oral_cavity_dose = dose_overrides.get(self.ModelRoiType.ORAL_CAVITY) if dose_overrides is not None and dose_overrides.get(self.ModelRoiType.ORAL_CAVITY) is not None else self.get_mean_roi_dose(self.ModelRoiType.ORAL_CAVITY)
        pcm_superior_dose = dose_overrides.get(self.ModelRoiType.PCM_SUPERIOR) if dose_overrides is not None and dose_overrides.get(self.ModelRoiType.PCM_SUPERIOR) is not None else self.get_mean_roi_dose(self.ModelRoiType.PCM_SUPERIOR)
        pcm_middle_dose = dose_overrides.get(self.ModelRoiType.PCM_MIDDLE) if dose_overrides is not None and dose_overrides.get(self.ModelRoiType.PCM_MIDDLE) is not None else self.get_mean_roi_dose(self.ModelRoiType.PCM_MIDDLE)
        pcm_inferior_dose = dose_overrides.get(self.ModelRoiType.PCM_INFERIOR) if dose_overrides is not None and dose_overrides.get(self.ModelRoiType.PCM_INFERIOR) is not None else self.get_mean_roi_dose(self.ModelRoiType.PCM_INFERIOR)

        exponent = {"oral_cavity":      constants["oral_cavity"][self._grade] * oral_cavity_dose,
                    "pcm_superior":     constants["pcm_superior"][self._grade] * pcm_superior_dose,
                    "pcm_middle":       constants["pcm_middle"][self._grade] * pcm_middle_dose,
                    "pcm_inferior":     constants["pcm_inferior"][self._grade] * pcm_inferior_dose}

        ntcp_total = 1 / (1 + math.exp(-intermediate_constant - sum(exponent.values())))

        return ntcp_total

# rs_ntcp_kno/viewmodels/mainviewmodel.py
class KNONTCP:

    def __init__(self, params) -> None:
        
        # Step 1 - Get the plans to compare to
        self._selected_plan1 = params.get(config.KEY_NTCP_PLAN1, None)
        self._selected_plan2 = params.get(config.KEY_NTCP_PLAN2, None)
        
        # Step 2 - Get RStation objects
        self._patient     = connect.get_current("Patient")
        self._case        = self._patient.Cases[0]
        if getRTPlanIndex(self._case, self._selected_plan1) > -1:
            self._plan    = self._case.TreatmentPlans[self._selected_plan1]
            self._roi_names = [roi.OfRoi.Name for roi in self._plan.GetTotalDoseStructureSet().RoiGeometries if roi.HasContours()]
            if "Glnds_Submand" not in self._roi_names:
                smdAlgebraStatus = self._do_smd_algebra()
                self._roi_names = [roi.OfRoi.Name for roi in self._plan.GetTotalDoseStructureSet().RoiGeometries if roi.HasContours()]

            self._roi_names_with_removed = [ModelRoi.REMOVED_NAME]
            self._roi_names_with_removed.extend(self._roi_names)
        else:
            print (' - [ERROR][KNONTCP] Plan 1 not found in plan list', self._selected_plan1)

        # Step 3 - Get model parameters (NOTE: this will be passed as params here)
        self._treatment_type        = params.get(config.KEY_NTCP_TREATMENT_TYPE, None)
        self._tumor_location        = params.get(config.KEY_NTCP_TUMOR_LOCATION, None)
        self._baseline_xerostomia   = params.get(config.KEY_NTCP_BASELINE_XEROSTOMIA, None)
        self._baseline_dysphagia    = params.get(config.KEY_NTCP_BASELINE_DYSPHAGIA, None)
        self._is_total_laryngectomy = params.get(config.KEY_NTCP_IS_TOTAL_LARYNGECTOMY, None)
        
        self._parotid_left   = "Parotid_L" if "Parotid_L" in self._roi_names else None
        self._parotid_right  = "Parotid_R" if "Parotid_R" in self._roi_names else None
        self._submandibulars = "Glnds_Submand" if "Glnds_Submand" in self._roi_names else None
        self._oral_cavity    = "Oral_Cavity" if "Oral_Cavity" in self._roi_names else None
        self._pcm_superior   = "Musc_Constrict_S" if "Musc_Constrict_S" in self._roi_names else None
        self._pcm_middle     = "Musc_Constrict_M" if "Musc_Constrict_M" in self._roi_names else None
        self._pcm_inferior   = "Musc_Constrict_I" if "Musc_Constrict_I" in self._roi_names else None
        self.ParotidsRemoved = params.get(config.KEY_PAROTIDS_REMOVED, False)
        if self._parotid_left is None:
            print (' - [INFO][KNONTCP] Parotid_L not found in ROI list')
        if self._parotid_right is None:
            print (' - [INFO][KNONTCP] Parotid_R not found in ROI list')
        if self._submandibulars is None:
            print (' - [INFO][KNONTCP] Glnds_Submand not found in ROI list')
        if self._oral_cavity is None:
            print (' - [INFO][KNONTCP] Oral_Cavity not found in ROI list')
        if self._pcm_superior is None:
            print (' - [INFO][KNONTCP] Musc_Constrict_S not found in ROI list')
        if self._pcm_middle is None:
            print (' - [INFO][KNONTCP] Musc_Constrict_M not found in ROI list')
        if self._pcm_inferior is None:
            print (' - [INFO][KNONTCP] Musc_Constrict_I not found in ROI list')

        # Step 4 - Define NTCP model objects
        plan1Index = getRTPlanIndex(self._case, self._selected_plan1)
        plan2Index = getRTPlanIndex(self._case, self._selected_plan2)

        if plan1Index > -1:
            if self._parotid_left is not None and self._parotid_right is not None and self._submandibulars is not None:
                self.plan1_ntcp_xerostomia_grade_2 = NtcpKnoXerostomia(self._case.TreatmentPlans[plan1Index], 2)
                self.plan1_ntcp_xerostomia_grade_3 = NtcpKnoXerostomia(self._case.TreatmentPlans[plan1Index], 3)
                self._set_xerostomia_variables(self.plan1_ntcp_xerostomia_grade_2)
                self._set_xerostomia_variables(self.plan1_ntcp_xerostomia_grade_3)
            else:
                self.plan1_ntcp_xerostomia_grade_2: NtcpKnoXerostomia = None
                self.plan1_ntcp_xerostomia_grade_3: NtcpKnoXerostomia = None    

            if self._oral_cavity is not None and self._pcm_superior is not None and self._pcm_middle is not None and self._pcm_inferior is not None:
                self.plan1_ntcp_dysphagia_grade_2  = NtcpKnoDysphagia(self._case.TreatmentPlans[plan1Index], 2)
                self.plan1_ntcp_dysphagia_grade_3  = NtcpKnoDysphagia(self._case.TreatmentPlans[plan1Index], 3)
                self._set_dysphagia_variables(self.plan1_ntcp_dysphagia_grade_2)
                self._set_dysphagia_variables(self.plan1_ntcp_dysphagia_grade_3)
            else:
                self.plan1_ntcp_dysphagia_grade_2 : NtcpKnoDysphagia  = None
                self.plan1_ntcp_dysphagia_grade_3 : NtcpKnoDysphagia  = None    
        else:
            self.plan1_ntcp_xerostomia_grade_2: NtcpKnoXerostomia = None
            self.plan1_ntcp_xerostomia_grade_3: NtcpKnoXerostomia = None
            self.plan1_ntcp_dysphagia_grade_2 : NtcpKnoDysphagia  = None
            self.plan1_ntcp_dysphagia_grade_3 : NtcpKnoDysphagia  = None
        
        if plan2Index > -1:
            self.plan2_ntcp_xerostomia_grade_2 = NtcpKnoXerostomia(self._case.TreatmentPlans[plan2Index], 2)
            self.plan2_ntcp_xerostomia_grade_3 = NtcpKnoXerostomia(self._case.TreatmentPlans[plan2Index], 3)
            self.plan2_ntcp_dysphagia_grade_2  = NtcpKnoDysphagia(self._case.TreatmentPlans[plan2Index], 2)
            self.plan2_ntcp_dysphagia_grade_3  = NtcpKnoDysphagia(self._case.TreatmentPlans[plan2Index], 3)
            self._set_xerostomia_variables(self.plan2_ntcp_xerostomia_grade_2)
            self._set_xerostomia_variables(self.plan2_ntcp_xerostomia_grade_3)
            self._set_dysphagia_variables(self.plan2_ntcp_dysphagia_grade_2)
            self._set_dysphagia_variables(self.plan2_ntcp_dysphagia_grade_3)
        else:
            self.plan2_ntcp_xerostomia_grade_2: NtcpKnoXerostomia = None
            self.plan2_ntcp_xerostomia_grade_3: NtcpKnoXerostomia = None
            self.plan2_ntcp_dysphagia_grade_2 : NtcpKnoDysphagia  = None
            self.plan2_ntcp_dysphagia_grade_3 : NtcpKnoDysphagia  = None

        # Step 99 - Results
        self._ntcp_xerostomia_grade_2_values = None
        self._ntcp_xerostomia_grade_3_values = None
        self._ntcp_dysphagia_grade_2_values = None
        self._ntcp_dysphagia_grade_3_values = None
        self._sum_grade_2_values = None
        self._sum_grade_3_values = None

    def _do_smd_algebra(self):
        
        import numpy as np
        algebraStatus = False

        def getMarginSettings(distInCm):
            return { 'Type': "Expand", 'Superior': distInCm, 'Inferior': distInCm, 'Anterior': distInCm, 'Posterior': distInCm, 'Right': distInCm, 'Left': distInCm }

        params = [
            {
                'roiNameNew': 'Glnds_Submand'
                , 'expARois': ['Glnd_Submand_L']
                , 'expAMarginSettings': getMarginSettings(0)
                , 'expBRois': ['Glnd_Submand_R']
                , 'expBMarginSettings': getMarginSettings(0)
                , 'roiType': 'Organ'
            }
        ]

        for param in params:
            try:
                roiNameNew = param['roiNameNew']
                expARois   = param['expARois']
                expBRois   = param['expBRois']
                roiNameNewBool = roiNameNew in self._roi_names
                expABool   = np.any(roiName in self._roi_names for roiName in expARois)
                expBBool   = np.any(roiName in self._roi_names for roiName in expBRois)
                if not roiNameNewBool:
                    if expABool and expBBool:
                        if roiNameNew not in self._roi_names:
                            print (f'\n  -- [doROIAlgebraForProtonAutoContours] Creating ROI: {roiNameNew} \n')
                            self._case.PatientModel.CreateRoi(Name=roiNameNew, Color="Green", Type=param['roiType'], TissueName=None, RbeCellTypeName=None, RoiMaterial=None)
                        self._case.PatientModel.RegionsOfInterest[roiNameNew].SetAlgebraExpression(
                            ExpressionA={ 'Operation': "Union", 'SourceRoiNames': expARois, 'MarginSettings': param['expAMarginSettings'] }
                            , ExpressionB={ 'Operation': "Union", 'SourceRoiNames': expBRois, 'MarginSettings': param['expBMarginSettings']}
                            , ResultOperation="Union", ResultMarginSettings=getMarginSettings(0)
                        )
                        self._case.PatientModel.RegionsOfInterest[roiNameNew].UpdateDerivedGeometry(Examination=self._case.Examinations[0], Algorithm="Auto")
                        roiNameNewObj = self._case.PatientModel.StructureSets[0].RoiGeometries[roiNameNew]
                        print ('')
                        if roiNameNewObj.HasContours():
                            print (f'  -- [doROIAlgebraForProtonAutoContours] ROI: {roiNameNew} created with volune: {roiNameNewObj.GetRoiVolume()}')
                            self._patient.Save()
                            algebraStatus = True
                        else:
                            print (f'  -- [doROIAlgebraForProtonAutoContours] ROI: {roiNameNew} created but has no contours')
                        print ('')

                    else:
                        print (f'  -- [doROIAlgebraForProtonAutoContours] expABool: {expABool} and expBBool:{expBBool}')
                else:
                    print (f'  -- [doROIAlgebraForProtonAutoContours] ROI: {roiNameNew} already exists')

            except:
                traceback.print_exc()
        
        return algebraStatus
    
    def _set_dysphagia_variables(self, model: NtcpKnoDysphagia):
        try:
            if model is None:
                return
            model.set_model_roi(self._oral_cavity, NtcpKnoDysphagia.ModelRoiType.ORAL_CAVITY)
            model.set_model_roi(self._pcm_superior, NtcpKnoDysphagia.ModelRoiType.PCM_SUPERIOR)
            model.set_model_roi(self._pcm_middle, NtcpKnoDysphagia.ModelRoiType.PCM_MIDDLE)
            model.set_model_roi(self._pcm_inferior, NtcpKnoDysphagia.ModelRoiType.PCM_INFERIOR)
            if self._baseline_dysphagia is not None:
                model.set_baseline(NtcpKnoDysphagia.BaselineType(self._baseline_dysphagia))
            if self._treatment_type is not None:
                model.set_treatment_type(TreatmentType(self._treatment_type))
            if self._tumor_location is not None:
                model.set_tumor_location(NtcpKnoDysphagia.TumorLocationType(self._tumor_location))

        except Exception as e:
            print(f" - [KNONTCP._set_dysphagia_variables()] Error in Dysphagia variables:\n{e}")

    def _set_xerostomia_variables(self, model: NtcpKnoXerostomia):
        try:
            if model is None:
                return
            model.set_model_roi(self._parotid_left, NtcpKnoXerostomia.ModelRoiType.PAROTID_LEFT)
            model.set_model_roi(self._parotid_right, NtcpKnoXerostomia.ModelRoiType.PAROTID_RIGHT)
            model.set_model_roi(self._submandibulars, NtcpKnoXerostomia.ModelRoiType.SUBMANDIBULARS)
            if self._baseline_xerostomia is not None:
                model.set_baseline(NtcpKnoXerostomia.BaselineType(self._baseline_xerostomia))
            if self._treatment_type is not None:
                model.set_treatment_type(TreatmentType(self._treatment_type))
        except Exception as e:
            print(f"Error in Xerostomia variables:\n{e}")

    def _get_ntcp_model(self, plan1_model: NtcpKnoAbstractBase, plan2_model: NtcpKnoAbstractBase):
        
        if plan1_model is None and plan2_model is None:
            return None
        values = {}
        
        try:
            values = {'ModelName':  plan1_model.model_name if plan1_model is not None else plan2_model.model_name,
                        'Plan_1':     None,
                        'Plan_2':     None,
                        'Criterion':  plan1_model.ntcp_proton_criterion if plan1_model is not None else plan2_model.ntcp_proton_criterion,
                        'NtcpPerGy':  None,
                        'Indication': None
                    }

            if plan1_model is not None:
                
                try:
                    plan1_model.plan.TreatmentCourse.TotalDose.UpdateDoseGridStructures() # done for dose.GetDoseStatistic(RoiName="Glnds_Submand", DoseType="Average") in get_mean_roi_dose
                except:
                    traceback.print_exc()

                values['Plan_1'] = {'Name':                plan1_model.plan.Name,
                                    'Modality':            plan1_model.plan.BeamSets[0].Modality,
                                    'PlanValue':           plan1_model.get_ntcp(),
                                    'BaseNtcp':            plan1_model.get_base_ntcp(),
                                    'Threshold':           plan1_model.get_ntcp_threshold(),
                                    'ComparisonJustified': None,
                                    'Rois':                {}
                                    }
                for roi_type, roi in plan1_model.model_rois.items():
                    values['Plan_1']['Rois'][roi.name] = {
                        'mean_dose': plan1_model.get_mean_roi_dose(roi_type), 'gradient': plan1_model.get_ntcp_per_gy(roi_type)
                    }
                if plan1_model.plan.BeamSets[0].Modality == "Photons":
                    values['Plan_1']['ComparisonJustified'] = values['Plan_1']['PlanValue'] > values['Plan_1']['Threshold']

            if plan2_model is not None:
                values['Plan_2'] = {'Name':           plan2_model.plan.Name,
                                    'Modality':       plan2_model.plan.BeamSets[0].Modality,
                                    'PlanValue':      plan2_model.get_ntcp(),
                                    'BaseNtcp':       plan2_model.get_base_ntcp(),
                                    'DeltaNtcp':      None,
                                    'Threshold':      None,
                                    'PlanWins':       False,
                                    'Rois':           {}
                                    }
                for roi_type, roi in plan2_model.model_rois.items():
                    values['Plan_2']['Rois'][roi.name] = {'mean_dose': plan2_model.get_mean_roi_dose(roi_type), 'gradient': plan2_model.get_ntcp_per_gy(roi_type), 'color': "Black"}

                if plan1_model is not None:
                    values['Plan_2']['DeltaNtcp'] = values['Plan_1']['PlanValue'] - values['Plan_2']['PlanValue']
                    values['Plan_2']['Threshold'] = plan2_model.get_ntcp_threshold()
                    values['Plan_2']['PlanWins'] = values['Plan_2']['DeltaNtcp'] >= plan1_model.ntcp_proton_criterion and values['Plan_1']['ComparisonJustified']
                    if plan1_model.plan.BeamSets[0].Modality == "Photons" and plan2_model.plan.BeamSets[0].Modality == "Protons":
                        values['Indication'] = plan2_model.plan.BeamSets[0].Modality if values['Plan_2']['PlanWins'] else plan1_model.plan.BeamSets[0].Modality

            return values
        except Exception as e:
            print(f"Error in creating NTCP values")
            traceback.print_exc()
            return values

    @property
    def ntcp_xerostomia_grade_2_values(self) -> dict:
        
        if self._parotid_left is None or self._parotid_right is None or self._submandibulars is None:
            return None
        
        if self.ParotidsRemoved:
            return None
        if self._ntcp_xerostomia_grade_2_values is None:
            self._ntcp_xerostomia_grade_2_values = self._get_ntcp_model(self.plan1_ntcp_xerostomia_grade_2, self.plan2_ntcp_xerostomia_grade_2)
        return self._ntcp_xerostomia_grade_2_values

    @property
    def ntcp_xerostomia_grade_3_values(self) -> dict:

        if self._parotid_left is None or self._parotid_right is None or self._submandibulars is None:
            return None
        
        if self.ParotidsRemoved:
            return None
        if self._ntcp_xerostomia_grade_3_values is None:
            self._ntcp_xerostomia_grade_3_values = self._get_ntcp_model(self.plan1_ntcp_xerostomia_grade_3, self.plan2_ntcp_xerostomia_grade_3)
        return self._ntcp_xerostomia_grade_3_values

    @property
    def ntcp_dysphagia_grade_2_values(self) -> dict:

        if self._oral_cavity is None or self._pcm_superior is None or self._pcm_middle is None or self._pcm_inferior is None:
            return None

        if self._is_total_laryngectomy:
            return None
        if self._ntcp_dysphagia_grade_2_values is None:
            self._ntcp_dysphagia_grade_2_values = self._get_ntcp_model(self.plan1_ntcp_dysphagia_grade_2, self.plan2_ntcp_dysphagia_grade_2)
        return self._ntcp_dysphagia_grade_2_values

    @property
    def ntcp_dysphagia_grade_3_values(self) -> dict:
        if self._is_total_laryngectomy:
            return None
        if self._ntcp_dysphagia_grade_3_values is None:
            self._ntcp_dysphagia_grade_3_values = self._get_ntcp_model(self.plan1_ntcp_dysphagia_grade_3, self.plan2_ntcp_dysphagia_grade_3)
        return self._ntcp_dysphagia_grade_3_values
