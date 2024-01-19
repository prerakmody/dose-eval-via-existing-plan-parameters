"""
This script expects
    - a folder containing CT, RTDose, RTPlan, RTStruct .dcms
    - assets/objective-template-photon-kno.xml
    - assets/eval-template-photon.csv
    - assets/isodose.xml

Run the script from "__main__"
"""

# Import RStation libraries
import connect

# Import private libraries
import helpers as helpers
import config as config

# Import general libraries
import re
import sys
import pdb
import json
import time
import shutil
import logging
import datetime
import traceback
import numpy as np
from pathlib import Path  

DEBUG_PDB = False
def print(*args, **kwargs):
    logging.info(" ".join(map(str, args)), **kwargs)
    

########################################################
#                        HELPERS                       #
########################################################

def optimizePlan(planName, count=10, reset=True, pathIsoDoseXML=None):

    optimizeStatus = False
    objectiveValues = []
    try:
        # Step 1 - Get plan and beamset
        patient = helpers.rayStationSave()
        case    = patient.Cases[0]
        plan    = case.TreatmentPlans[planName]
        beamset = case.TreatmentPlans[planName].BeamSets[planName]
        beamSetIndex = plan.BeamSets.IndexOf(beamset)

        # Step 2 - Reset existing optimization
        if reset:
            plan.PlanOptimizations[beamSetIndex].ResetOptimization()
            print ('\n\n - [optimizePlan()][Patient={}][Plan={}] Optimization has been reset ...'.format(helpers.getPatientIdentifier(patient), planName))
        else:
            print ('\n\n - [optimizePlan()][Patient={}][Plan={}] Optimization has not been reset ...'.format(helpers.getPatientIdentifier(patient), planName))
            if plan.PlanOptimizations[beamSetIndex].ProgressOfOptimization is not None:
                objectiveValue = plan.PlanOptimizations[beamSetIndex].ProgressOfOptimization.ObjectiveValues[-1] # [TODO: range = (0.36,?)]
                print (' - [optimizePlan()][Patient={}][Plan={}] Previous objective {:.4f}'.format(helpers.getPatientIdentifier(patient), planName, objectiveValue))
        _ = helpers.rayStationSave()

        # Step 3 - Run optimization
        times = []
        
        for runID in range(count):
            t0 = time.time()
            print (' - [optimizePlan()][Patient={}][Plan={}] Running optimization step {}/{} ... '.format(helpers.getPatientIdentifier(patient), planName, runID+1, count))
            plan.PlanOptimizations[beamSetIndex].RunOptimization()   
            times.append(time.time() - t0)
            try:
                objectiveValue = plan.PlanOptimizations[beamSetIndex].ProgressOfOptimization.ObjectiveValues[-1] # [TODO: range = (0.36,?)]
                print (' --- [optimizePlan()] Optimization step {}/{} took {:.2f} seconds with mean objective: {:.4f}'.format(runID+1, count, times[-1], objectiveValue))
                objectiveValues.append(objectiveValue)
            except:
                traceback.print_exc()
                if DEBUG_PDB: pdb.set_trace()
        print (' - [optimizePlan()] Optimization took total {:.2f} seconds with mean objective: {:.4f}'.format(np.sum(times), objectiveValue))

        if runID == count-1:
            optimizeStatus = True
        
        if optimizeStatus:
            helpers.applyIsoDoseColors(pathIsoDoseXML, case, beamset)
        else:
            print (' - [optimizePlan()] Full optimization failed, not applying isodose colors')

    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace() 
    
    return optimizeStatus, objectiveValues

def uploadORUpdateObjectives(planName, pathKNOObjectives, uploadObjectivesBool, updateObjectivesBool, forceObjectives, objectiveFType):

    # Step 1 - Init
    objectiveStatus    = False

    # Step 2 - Boilerplate code for get patient and plans
    _, case, plan, beamset  = helpers.getPatientAndPlan(planName)

    # Step 3 - Get objectives
    objectivesFromPath = helpers.getObjectivesFromPath(plan, beamset, pathKNOObjectives)
    objectivesFromRS   = helpers.getObjectivesFromPlan(plan, beamset)

    # Step 4 - Upload/Update Objectives
    if uploadObjectivesBool:
        if len(objectivesFromRS) == 0 or forceObjectives:
            objectiveStatus = helpers.uploadObjectivesToRS(plan, beamset, objectivesFromPath)
        else:
            print(' - [uploadORUpdateObjectives()] Objectives already exist, not uploading')
    
    elif updateObjectivesBool:
        if len(objectivesFromRS) > 0:
            objectiveStatus = updateObjectives(planName, objectivesFromRS, objectivesFromPath, objectiveFType)
        else:
            print(' - [uploadORUpdateObjectives()] No objectives found, not updating')
    
    # Step 5 - Do sanity check
    helpers.doPlanSanityCheck(case, planName)

    return objectiveStatus

def checkIfRoiIgnore(planName, objectiveFType, roiName):
    
    # Step 1 - Init
    roisFullIgnore = []
    roisPartialIgnore = []
    roiIgnoreBool = False

    keyPlanCSTerm        = config.SUFFIX_PLAN_CS.format(config.PREFIX_CLINICAL_CONTOURS)
    keyPlanCSAutoTerm    = config.SUFFIX_PLAN_CS.format(config.PREFIX_AUTOMATED_CONTOURS)
    keyPlanDFOTerm       = config.SUFFIX_PLAN_DFO.format(config.PREFIX_CLINICAL_CONTOURS)
    keyPlanDFOAutoTerm   = config.SUFFIX_PLAN_DFO.format(config.PREFIX_AUTOMATED_CONTOURS)
    keyPlanDFO2Term      = config.SUFFIX_PLAN_DFO2.format(config.PREFIX_CLINICAL_CONTOURS)
    keyPlanDFO2AutoTerm  = config.SUFFIX_PLAN_DFO2.format(config.PREFIX_AUTOMATED_CONTOURS)
    keyPlanEUDTerm       = config.SUFFIX_PLAN_EUD.format(config.PREFIX_CLINICAL_CONTOURS)
    keyPlanEUDAutoTerm   = config.SUFFIX_PLAN_EUD.format(config.PREFIX_AUTOMATED_CONTOURS)
    keyPlanFinalTerm     = config.SUFFIX_PLAN_FINAL.format(config.PREFIX_CLINICAL_CONTOURS)
    keyPlanFinalAutoTerm = config.SUFFIX_PLAN_FINAL.format(config.PREFIX_AUTOMATED_CONTOURS)

    # Step 2 - Ignore check on certain ROIs
    if roiName in [config.KEYNAME_BRAINSTEM, config.KEYNAME_BRAINSTEM + config.KEY_3MM_SUFFIX
                   , config.KEYNAME_SPINALCORD, config.KEYNAME_SPINALCORD + config.KEY_3MM_SUFFIX
                   , config.KEYNAME_GHOST_CRANIAL, config.KEYNAME_EAR_L_GHOST, config.KEYNAME_EAR_R_GHOST
                   ]:
        return roiIgnoreBool

    # Step 2 - Get ROIs to ignore (while updating objectives)
    roisPartialIgnoreList = [config.REGEX_PROSTHESE, config.REGEX_DFO, config.REGEX_DMAX, config.REGEX_DOSE, config.REGEX_KAAK
                             , config.REGEX_LIPPEN
                             , config.KEYNAME_MANDIBLE_PTV
                             , config.KEYNAME_MANDIBLE
                             , config.REGEX_5805
                             , config.KEYNAME_BRAIN
                            ]
    
    if (keyPlanCSTerm in planName or keyPlanCSAutoTerm in planName) and objectiveFType == None: # this will actually never be called
        roisFullIgnore = [config.KEYNAME_BODY, config.KEYNAME_RING_PTV_DL2, config.KEYNAME_GHOST]

    elif (keyPlanDFOTerm in planName or keyPlanDFOAutoTerm in planName) and objectiveFType == config.KEY_FTYPE_DOSEFALLOFF:
        roisFullIgnore = [config.KEYNAME_BODY, config.KEYNAME_RING_PTV_DL2, config.KEYNAME_GHOST]
        roisPartialIgnore = roisPartialIgnoreList

    elif (keyPlanDFO2Term in planName or keyPlanDFO2AutoTerm in planName) and objectiveFType == config.KEY_FTYPE_DOSEFALLOFF:
        roisFullIgnore = [config.KEYNAME_BODY, config.KEYNAME_RING_PTV_DL2, config.KEYNAME_GHOST]
        roisPartialIgnore = roisPartialIgnoreList
    
    elif (keyPlanEUDTerm in planName or keyPlanEUDAutoTerm in planName) and objectiveFType == config.KEY_FTYPE_MAXEUD:
        roisFullIgnore = [config.KEYNAME_BODY, config.KEYNAME_RING_PTV_DL2, config.KEYNAME_GHOST]
        roisPartialIgnore = roisPartialIgnoreList
    
    elif (keyPlanFinalTerm in planName or keyPlanFinalAutoTerm in planName) and objectiveFType == None:
        roisFullIgnore = [config.KEYNAME_BODY]
        roisPartialIgnore = []

    # Step 3 - Check if ROI should be ignored
    for roiFullIgnore in roisFullIgnore:
        if roiName == roiFullIgnore:
            roiIgnoreBool = True
            break
    
    for roiPartialIgnore in roisPartialIgnore:
        if re.search(roiPartialIgnore, str(roiName).lower()) is not None:
            roiIgnoreBool = True
            break

    return roiIgnoreBool

def updateObjectiveRingPTVDL1(planName, objectivesFromRS, rsIdx):

    # Step 0 - Init
    updateStatus = False 
    existingRoiName      = objectivesFromRS[rsIdx].ForRegionOfInterest.Name
    existingWeight       = objectivesFromRS[rsIdx].DoseFunctionParameters.Weight

    keyPlanDFOTerm       = config.SUFFIX_PLAN_DFO.format(config.PREFIX_CLINICAL_CONTOURS)
    keyPlanDFOAutoTerm   = config.SUFFIX_PLAN_DFO.format(config.PREFIX_AUTOMATED_CONTOURS)
    keyPlanDFO2Term      = config.SUFFIX_PLAN_DFO2.format(config.PREFIX_CLINICAL_CONTOURS)
    keyPlanDFO2AutoTerm  = config.SUFFIX_PLAN_DFO2.format(config.PREFIX_AUTOMATED_CONTOURS)
    
    if (keyPlanDFOTerm in planName or keyPlanDFOAutoTerm in planName):
        # newWeight = 40 # [40, 50]
        # objectivesFromRS[rsIdx].DoseFunctionParameters.Weight = newWeight
        updateStatus = True
        print (f' - [updateObjectiveRingPTVDL1()] \troi: {existingRoiName}, \tweight: {existingWeight} --> {objectivesFromRS[rsIdx].DoseFunctionParameters.Weight}')
        print (f' - [updateObjectiveRingPTVDL1()] Not updating! ')

    elif (keyPlanDFO2Term in planName or keyPlanDFO2AutoTerm in planName):
        # newWeight = 80 # [80, 100]  
        # objectivesFromRS[rsIdx].DoseFunctionParameters.Weight = newWeight
        updateStatus = True
        print (f' - [updateObjectiveRingPTVDL1()] \troi: {existingRoiName}, \tweight: {existingWeight} --> {objectivesFromRS[rsIdx].DoseFunctionParameters.Weight}')
        print (f' - [updateObjectiveRingPTVDL1()] Not updating! ')
        
    else:
        pass

    return updateStatus
    
def updateObjectives(planName, objectivesFromRS, objectivesFromPath, objectiveFType):
    """
    Called from uploadORUpdateObjectives()

    Params:
    -------
    planName: str
    objectivesFromRS: list
    objectivesFromPath: list
    objectiveFType: str, ['DoseFallOff', 'DoseFallOff', 'MaxEUD', None]
    """
    updateObjectiveStatus = False

    try:

        # Step 1 - Init
        patient, _, _, _  = helpers.getPatientAndPlan(planName)
        keyPlanDFO2Term      = config.SUFFIX_PLAN_DFO2.format(config.PREFIX_CLINICAL_CONTOURS)
        keyPlanDFO2AutoTerm  = config.SUFFIX_PLAN_DFO2.format(config.PREFIX_AUTOMATED_CONTOURS)
        keyPlanFinalTerm     = config.SUFFIX_PLAN_FINAL.format(config.PREFIX_CLINICAL_CONTOURS) # R5 
        keyPlanFinalAutoTerm = config.SUFFIX_PLAN_FINAL.format(config.PREFIX_AUTOMATED_CONTOURS) # A5

        # Step 2 - Get existing objectives
        updatedFromPathObjectivesStatus = {}
        for pathIdx, objectiveFromPath in enumerate(objectivesFromPath):
            newRoiName      = objectiveFromPath.roi_name
            if newRoiName == config.KEYNAME_BODY: continue
            newFunctionType = objectiveFromPath.function_type
            key = f'{newRoiName}-{newFunctionType}'
            updatedFromPathObjectivesStatus[key] = {'exists': False, 'idx':pathIdx}

        # Step 4 - Update objectives
        if len(objectivesFromRS) > 0:
            print ('\n\n - [INFO][updateObjectives()][Patient={}] Updating Objectives ... '.format(helpers.getPatientIdentifier(patient)))
            
            # Step 4.1 - Update weights of existing objectives (based on function-type and roi-name)
            for rsIdx, objectiveFromRS in enumerate(objectivesFromRS):
                try:
                    existingRoiName      = objectiveFromRS.ForRegionOfInterest.Name
                    existingWeight       = objectiveFromRS.DoseFunctionParameters.Weight
                    if existingRoiName == config.KEYNAME_BODY: continue
                    try:
                        existingFunctionType = objectiveFromRS.DoseFunctionParameters.FunctionType
                    except:
                        existingFunctionType = config.KEY_FTYPE_DOSEFALLOFF
                    
                    # print (' - existingRoiName: {}/{}, existingFunctionType: {}, existingWeight: {}'.format(existingRoiName, config.KEYNAME_RING_LT_PTV_DL1, existingFunctionType, existingWeight))
                    if existingRoiName == config.KEYNAME_RING_LT_PTV_DL1 and objectiveFType is not None:
                        if updateObjectiveRingPTVDL1(planName, objectivesFromRS, rsIdx):
                            continue
                    
                    # [DEBUG]
                    if 0:
                        print (' - [DEBUG][updateObjectives()] roi: {}, \tfunctionType: {}, \tweight: {}'.format(existingRoiName, existingFunctionType, existingWeight))
                    
                    # Make note of existing objectives that are in new objectives
                    key = f'{existingRoiName}-{existingFunctionType}'
                    if key in updatedFromPathObjectivesStatus: updatedFromPathObjectivesStatus[key]['exists'] = True

                    # Step 4.1 - Check if existing objective is in new objectives
                    for objectiveFromPath in objectivesFromPath:
                        newRoiName      = objectiveFromPath.roi_name
                        newFunctionType = objectiveFromPath.function_type
                        newWeight       = objectiveFromPath.weight
                        if newRoiName == config.KEYNAME_BODY: continue
                        functionMatchBool = newFunctionType == objectiveFType
                        if  keyPlanFinalTerm in planName or keyPlanFinalAutoTerm in planName:
                            functionMatchBool = True # i.e update all rois

                        if functionMatchBool:
                            if not checkIfRoiIgnore(planName, objectiveFType, newRoiName):
                                if existingRoiName == newRoiName and existingFunctionType == newFunctionType:

                                    # Step 4.1.1 - Update weight 
                                    objectivesFromRS[rsIdx].DoseFunctionParameters.Weight = newWeight
                                    print (f'  -- [updateObjectives()] \troi: {existingRoiName}, \tfType: {existingFunctionType}, \tweight: {existingWeight} --> {newWeight}')

                                    # Step 4.1.2 - Update dose level(s)
                                    if existingFunctionType == config.KEY_FTYPE_MAXEUD:
                                        existingDoseLevel = objectivesFromRS[rsIdx].DoseFunctionParameters.DoseLevel
                                        existingEudParameterA = objectivesFromRS[rsIdx].DoseFunctionParameters.EudParameterA
                                        objectivesFromRS[rsIdx].DoseFunctionParameters.DoseLevel     = objectiveFromPath.doselevel
                                        objectivesFromRS[rsIdx].DoseFunctionParameters.EudParameterA = objectiveFromPath.eud_parameter_a
                                        print (f' --- [updateObjectives()] \tDoseLevel: {existingDoseLevel} --> {objectiveFromPath.doselevel}')
                                        print (f' --- [updateObjectives()] \tEudParameterA: {existingEudParameterA} --> {objectiveFromPath.eud_parameter_a}')
                                    elif existingFunctionType == config.KEY_FTYPE_DOSEFALLOFF:
                                        if keyPlanDFO2Term in planName or keyPlanDFO2AutoTerm in planName: # [Ref: https://iprova.lumc.nl/Portal/#/document/1fc74366-40c1-4b84-a653-0455b6c891f8 (Vervolgens voor beide opties)]
                                            print (f' --- [updateObjectives()] \tHighDoseLevel: {objectiveFromRS.DoseFunctionParameters.HighDoseLevel} --> {objectiveFromPath.high_doselevel}')
                                            print (f' --- [updateObjectives()] \tLowDoseLevel : {objectiveFromRS.DoseFunctionParameters.LowDoseLevel}  --> {objectiveFromPath.low_doselevel}')
                                            objectivesFromRS[rsIdx].DoseFunctionParameters.HighDoseLevel = objectiveFromPath.high_doselevel
                                            objectivesFromRS[rsIdx].DoseFunctionParameters.LowDoseLevel  = objectiveFromPath.low_doselevel
                                
                                key = f'{newRoiName}-{newFunctionType}'
                                updatedFromPathObjectivesStatus[key]['updated'] = True
                except:
                    traceback.print_exc()
                    if DEBUG_PDB: pdb.set_trace()

            # Step 4.2 - Add new objectives (if 1) they dont exist, and if they match 2) function-type and 3) roiName conditions)
            print ('')
            for key in updatedFromPathObjectivesStatus:
                if updatedFromPathObjectivesStatus[key]['exists'] is False:
                    roiName = key.split('-')[0]
                    fType   = key.split('-')[1]
                    
                    fTypeCondition     = fType == objectiveFType # only if ftType in [DoseFallOff, MaxEUD]
                    roiIgnoreCondition = checkIfRoiIgnore(planName, objectiveFType, roiName) # this roiName comes from existing objectives of RS
                    if objectiveFType == None: # since this is the last plan, we want to add all (non-existent) objectives
                        fTypeCondition = True

                    if fTypeCondition and not roiIgnoreCondition:
                        print (f' - [updateKNOObjectivesinRStation()] Adding new objective for {roiName} of type {fType}')
                        newObjective = objectivesFromPath[updatedFromPathObjectivesStatus[key]['idx']]
                        print (f' --- [updateKNOObjectivesinRStation()] newObjective for roi: {newObjective.roi_name} and fType: {newObjective.function_type} with weight: {newObjective.weight}')
                        if newObjective.weight > 0:
                            newObjective.apply()
        
        else:
            print (f' - [updateKNOObjectivesinRStation()] No objectives found in RS for plan {planName}')
        
        updateObjectiveStatus = True

    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace()
    
    # DEBUG (for config.KEYNAME_RING_LT_PTV_DL1, ghost)
    if 1:
        for rsIdx, objectiveFromRS in enumerate(objectivesFromRS):
            try:
                existingRoiName      = str(objectiveFromRS.ForRegionOfInterest.Name)
                existingWeight       = objectiveFromRS.DoseFunctionParameters.Weight
                try:
                    existingFunctionType = objectiveFromRS.DoseFunctionParameters.FunctionType
                except:
                    existingFunctionType = config.KEY_FTYPE_DOSEFALLOFF
                if existingRoiName == config.KEYNAME_RING_LT_PTV_DL1:
                    print ('')
                    print (' - [DEBUG2][plan={}] existingRoiName: {}, existingFunctionType: {}, existingWeight: {}'.format(planName, existingRoiName, existingFunctionType, existingWeight))
                    print ('')
                
                elif config.KEYNAME_GHOST in existingRoiName.lower() or config.REGEX_DFO in existingRoiName.lower():
                    if existingRoiName not in [config.KEYNAME_GHOST_CRANIAL, config.KEYNAME_EAR_L_GHOST, config.KEYNAME_EAR_R_GHOST]:
                        print ('')
                        print (' - [DEBUG2][plan={}] existingRoiName: {}, existingFunctionType: {}, existingWeight: {}'.format(planName, existingRoiName, existingFunctionType, existingWeight))
                        print ('')
            except:
                pass


    return updateObjectiveStatus

def airOverride(setMaterial=True):

    try:
        case = connect.get_current(config.KEYNAME_CASE)
        REGEX_AIR_ROI = 'override lucht'

        # Step 1 - Check if override lucht exists
        roiExists = False
        roiName = None
        for roi in case.PatientModel.RegionsOfInterest:
            # print (roi.Name, REGEX_AIR_ROI)
            if re.search(REGEX_AIR_ROI, roi.Name.lower()) is not None:
                roiExists = True
                roiName = roi.Name
                break
                
        # Step 2 - If yes, then override material
        print ('\n\n ======================================================================== \n\n')
        print (' - [airOverride()] roiExists: ', roiExists)
        if roiExists and setMaterial:
            case.PatientModel.RegionsOfInterest[roiName].SetRoiMaterial(Material=case.PatientModel.Materials[1])
            sys.stdout.write(f' - [airOverride()] Overriding {roiName} with material: {case.PatientModel.Materials[1].Name}')
            print (f' - [airOverride()] Overriding {roiName} with material: {case.PatientModel.Materials[1].Name}')
        else:
            print (' - [airOverride()] No override lucht found \n\n')
        print ('\n\n ======================================================================== \n\n')

    except:
        traceback.print_exc()

########################################################
#                          NTCP                        #
########################################################

def getNTCPVals(patientID, plans, pathPatient):

    try:
        
        print (f' - [getNTCPVals()] Getting NTCP values for patient: {patientID}')
        
        # Step 1 - Other params
        params = {
            config.KEY_NTCP_PLAN1  : None
            , config.KEY_NTCP_PLAN2: None
            , config.KEY_NTCP_TREATMENT_TYPE     : helpers.TreatmentType.PRIMARY
            , config.KEY_NTCP_TUMOR_LOCATION     : helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX # NtcpKnoDysphagia.TumorLocationType.{ORAL_CAVITY, PHARYNX, LARYNX}
            , config.KEY_NTCP_BASELINE_XEROSTOMIA: helpers.NtcpKnoXerostomia.BaselineType.NONE # NtcpKnoXerostomia.BaselineType.{NONE, LITTLE, SEVERE}
            , config.KEY_NTCP_BASELINE_DYSPHAGIA : helpers.NtcpKnoDysphagia.BaselineType.GRADE_0_1 # NtcpKnoDysphagia.BaselineType.{GRADE_0_1, GRADE_2, GRADE_3_4}
            , config.KEY_NTCP_IS_TOTAL_LARYNGECTOMY: False
            , config.KEY_PAROTIDS_REMOVED          : False
        }

        # Step 2 - main
        res = {}
        if params[config.KEY_NTCP_TUMOR_LOCATION] is not None:
            for plan in plans:

                # Step 2.1 - Init
                params[config.KEY_NTCP_PLAN1] = plan
                res[plan] = {}

                initStatus = False
                try:
                    objNTCP = helpers.KNONTCP(params)
                    initStatus = True
                except:
                    traceback.print_exc()
                    res[plan] = {
                        config.KEY_XERO_GRADE2: -1
                        , config.KEY_XERO_GRADE3: -1
                        , config.KEY_DYS_GRADE2 : -1
                        , config.KEY_DYS_GRADE3 : -1
                    }

                # Step 2.2 - Calculate
                if initStatus:
                    try:
                        tmpX2 = objNTCP.ntcp_xerostomia_grade_2_values
                        if tmpX2 is not None:
                            if tmpX2[config.KEY_NTCP_PLAN_1] is not None:
                                res[plan][config.KEY_XERO_GRADE2] = round(tmpX2[config.KEY_NTCP_PLAN_1][config.KEY_NTCP_PLAN_VALUE], 5)
                            else:
                                res[plan][config.KEY_XERO_GRADE2] = -1    
                        else:
                            res[plan][config.KEY_XERO_GRADE2] = -1
                    except:
                        res[plan][config.KEY_XERO_GRADE2] = -1 

                    try:
                        tmpX3 = objNTCP.ntcp_xerostomia_grade_3_values
                        if tmpX3 is not None:
                            if tmpX3[config.KEY_NTCP_PLAN_1] is not None:
                                res[plan][config.KEY_XERO_GRADE3] = round(tmpX3[config.KEY_NTCP_PLAN_1][config.KEY_NTCP_PLAN_VALUE], 5)
                            else:
                                res[plan][config.KEY_XERO_GRADE3] = -1    
                        else:
                            res[plan][config.KEY_XERO_GRADE3] = -1
                    except:
                        res[plan][config.KEY_XERO_GRADE3] = -1
                    
                    try:
                        tmpD2 = objNTCP.ntcp_dysphagia_grade_2_values
                        if tmpD2 is not None:
                            if tmpD2[config.KEY_NTCP_PLAN_1] is not None:
                                res[plan][config.KEY_DYS_GRADE2] = round(tmpD2[config.KEY_NTCP_PLAN_1][config.KEY_NTCP_PLAN_VALUE], 5)
                            else:
                                res[plan][config.KEY_DYS_GRADE2] = -1                    
                        else:
                            res[plan][config.KEY_DYS_GRADE2] = -1
                    except:
                        res[plan][config.KEY_DYS_GRADE2] = -1

                    try:
                        tmpD3 = objNTCP.ntcp_dysphagia_grade_3_values
                        if tmpD3 is not None:
                            if tmpD3[config.KEY_NTCP_PLAN_1] is not None:
                                res[plan][config.KEY_DYS_GRADE3] = round(tmpD3[config.KEY_NTCP_PLAN_1][config.KEY_NTCP_PLAN_VALUE], 5)
                            else:
                                res[plan][config.KEY_DYS_GRADE3] = -1
                        else:
                            res[plan][config.KEY_DYS_GRADE3] = -1
                    except:
                        res[plan][config.KEY_DYS_GRADE3] = -1
                
        # Step 3 - Save
        # print (f' - [getNTCPVals()] res: {res}')
        pathPatientNTCP = Path(pathPatient) / config.FILENAME_NTCP_PHOTON_RESULTS
        with open(str(pathPatientNTCP), 'w', encoding='utf-8') as fp:
            json.dump(res, fp, indent=4)
        
        print (f' - [getNTCPVals()] Done NTCP calculations for {patientID}')
        
    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace()

########################################################
#                    AUTO-HELPERS                      #
########################################################
def doAutoContouring():

    autoContourStatus = False
    t0 = time.time()
    try:
        
        # Step 0 - Init
        print (f' \n\n ===================== start autocontouring ===================== \n\n')
        patient = helpers.rayStationSave()
        case    = patient.Cases[0]

        # Step 1.1 - Get OARs to include    
        oarDuplicateStatus = helpers.checkOARDuplicateStatus(case)
        oarsToInclude = []
        for oar in oarDuplicateStatus:
            if not oarDuplicateStatus[oar]:
                oarsToInclude.append(oar)
        print (' - [doAutoContouring()] OARs pending to auto-contour: ', oarsToInclude)
                
        # Step 1.2 - Run OAR segmentation
        if len(oarsToInclude):
            _ = helpers.rayStationSave()
            case.SetCurrent() # Why == Potential Error: "The case to which the examination belongs mist be selected"
            examination = case.Examinations[0]
            _ = examination.RunOarSegmentation(ModelName="RSL Head and Neck CT", ExaminationsAndRegistrations={ 'CT 1': None }, RoisToInclude=oarsToInclude)
            helpers.rayStationSave()
        else:
            print (' - [doAutoContouring()] No OARs to auto-contour')
                
        # Step 1.3 - Check aut--contours status
        ## If they dont have the config.KEY_AUTOCONTOUR_SUFFIX, then they were not initially present, so simply rename them
        oarDuplicateStatus = helpers.checkOARDuplicateStatus(case)
        for oar in oarDuplicateStatus:
            if not oarDuplicateStatus[oar]:
                print (' - [doAutoContouring()] OAR auto-contour not present: ', oar)
                if oar == config.KEYNAME_CAVITY_ORAL:
                    newName                                            = config.KEYNAME_ORAL_CAVITY + config.KEY_AUTOCONTOUR_SUFFIX
                    case.PatientModel.RegionsOfInterest[oar].Name      = newName
                    case.PatientModel.RegionsOfInterest[newName].Color = config.OAR_DUPLICATE_COLOR_RGB_STRING 
                    print (' - [doAutoContouring()] Renamed OAR: ', oar, ' --> ', newName)
                elif oar == config.KEYNAME_ESOPHAGUS_S:
                    try:
                        newName                                            = config.KEYNAME_ESOPHAGUS + config.KEY_AUTOCONTOUR_SUFFIX
                        case.PatientModel.RegionsOfInterest[oar].Name      = newName
                        case.PatientModel.RegionsOfInterest[newName].Color = config.OAR_DUPLICATE_COLOR_RGB_STRING 
                        print (' - [doAutoContouring()] Renamed OAR: ', oar, ' --> ', newName)
                    except:
                        traceback.print_exc()
                else:
                    case.PatientModel.RegionsOfInterest[oar].Name = oar + config.KEY_AUTOCONTOUR_SUFFIX

        # Step 1.4 - Do ROI Algebra
        helpers.doROIAlgebraForAutoContours(case)
        _ = helpers.rayStationSave()
        timeTaken = round(time.time() - t0, 2)
        print (f' \n\n ===================== end autocontouring (in {timeTaken} s) ===================== \n\n')
        return True, timeTaken

    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace()

    timeTaken = round(time.time() - t0, 2)
    return False, timeTaken

########################################################
#                        MAIN(S)                       #
########################################################

# Func 1
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
    helpers.rayStationSave()
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
                    else          : patient = helpers.getPatientById(patientIdCheck, lastFind=True)
                    if patient is None:
                        print ('\n - [uploadRTAppsDataToRStation()] --------------------- Step 1.1: Uploading patient data for patientId={} \n\n'.format(patientIdCheck))
                        if Path(pathPatientCTFolder).exists():
                            patientID, studyUID, seriesUID = helpers.updateCTDicoms(pathPatientCTFolder)
                            print ('\n - [uploadRTAppsDataToRStation()] --------------------- Step 1.1: Uploading CT data for {} ... \n\n'.format(patientID))
                            warnings = db.ImportPatientFromPath(Path=str(pathPatientCTFolder), SeriesOrInstances=[{'PatientID': patientID, 'StudyInstanceUID': str(studyUID), 'SeriesInstanceUID': str(seriesUID)}], ImportFilter='', BrachyPlanImportOverrides={})
                            patient = helpers.rayStationSave()
                            # patient.Cases[0].SetCurrent()
                            print ('\n - [uploadRTAppsDataToRStation()] --------------------- Step 1.1: Uploaded CT data for {} \n\n'.format(helpers.getPatientIdentifier(patient)))
                            helpers.setEquipmentName()
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
            patient = helpers.rayStationSave()
            if patient is not None and patientCTBool:
                case      = patient.Cases[0]
                casename  = case.CaseName
                patientID = patient.PatientID
                
                # Step 1.2 - Upload (existing) RTStruct
                try:
                    print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: Upload RTStruct data \n\n'.format(helpers.getPatientIdentifier(patient)))
                    if not helpers.checkForRTStruct(case):
                        pathPatientRTStructFile = [each for each in pathPatientRTStructFolder.iterdir()][0]
                        if Path(pathPatientRTStructFile).exists():
                            studyUIDRTStruct, seriesUIDRTStruct = helpers.updateRTStructDicoms(pathPatientRTStructFile)
                            print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: Uploading RTStruct data ... \n\n'.format(helpers.getPatientIdentifier(patient)))
                            warningsRTStruct  = patient.ImportDataFromPath(Path=str(pathPatientRTStructFolder), CaseName=casename, SeriesOrInstances=[{'PatientID': patientID, 'StudyInstanceUID': str(studyUIDRTStruct), 'SeriesInstanceUID': str(seriesUIDRTStruct)}], ImportFilter='', BrachyPlanImportOverrides={}, AllowMismatchingPatientID=True)
                            print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: Uploaded RTStruct data \n\n'.format(helpers.getPatientIdentifier(patient)))
                            patientRTStructBool = True
                            helpers.rayStationSave()
                        else:
                            print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: No RTStruct data found: {} \n\n'.format(helpers.getPatientIdentifier(patient), pathPatientRTStructFile))
                    else:
                        print ('\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: RTStruct data already exists \n\n'.format(helpers.getPatientIdentifier(patient)))
                        patientRTStructBool = True
                except:
                    traceback.print_exc()
                    print ('\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.2: Issue with RTStruct data upload: {} \n\n'.format(helpers.getPatientIdentifier(patient), pathPatientRTStructFile))
                    patientRTStructBool = False

                print ('\n ---------------------------------------------------------- ')
                    
                # Step 1.3 - Upload (existing) RTPlan/RTDose
                if patientCTBool and patientRTStructBool:
                    if not helpers.checkForRTPlan(case, planName):
                        print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.3: Upload RTPlan/RTDose data \n\n'.format(helpers.getPatientIdentifier(patient)))
                        studyUIDRTPlan  = Path(pathPatientRTPlanFolder).parts[-2]
                        seriesUIDRTPlan = Path(pathPatientRTPlanFolder).parts[-1].split('_')[-1]
                        studyUIDRTDose  = Path(pathPatientRTDoseFolder).parts[-2]
                        seriesUIDRTDose = Path(pathPatientRTDoseFolder).parts[-1].split('_')[-1]
                        
                        studyID = helpers.updateRTPlanDicoms(pathPatientRTPlanFolder)
                        # studyUIDRTPlan = studyID 
                        # studyUIDRTDose = studyID
                        
                        pathTempRTDoseAndRTPlanFolder = helpers.getTempRTDoseAndRTPlanFolder(pathPatientRTPlanFolder, pathPatientRTDoseFolder)
                        print ('\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.3: Uploading RTPlan/RTDose data ... \n\n'.format(helpers.getPatientIdentifier(patient)))
                        warningsRTPlan  = patient.ImportDataFromPath(Path=str(pathTempRTDoseAndRTPlanFolder), CaseName=casename, SeriesOrInstances=[
                            {'PatientID': patientID, 'StudyInstanceUID': str(studyUIDRTPlan), 'SeriesInstanceUID': str(seriesUIDRTPlan)}    
                            , {'PatientID': patientID, 'StudyInstanceUID': str(studyUIDRTDose), 'SeriesInstanceUID': str(seriesUIDRTDose)}
                            ], ImportFilter='', BrachyPlanImportOverrides={}, AllowMismatchingPatientID=True)
                        shutil.rmtree(pathTempRTDoseAndRTPlanFolder)
                        patientRTPlanBool = True
                        patientRTDoseBool = True
                        print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.3: Uploaded RTPlan/RTDose data \n\n'.format(helpers.getPatientIdentifier(patient)))
                    else:
                        print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.3: RTPlan/RTDose already exists \n\n'.format(helpers.getPatientIdentifier(patient)))
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

# Func 2
def copyPlanAndOptimize(basePlanName, newPlanName
                        , pathKNOObjectives, uploadObjectivesBool, updateObjectivesBool, forceObjectives, objectiveFType
                        , optSteps, optReset, pathIsoDoseXML=None
                        , debug=False):

    t0 = time.time()
    optimizeValue = -1
    timeTaken = -1
    try:

        # Step 0 - for console printing (while debugging across multiple consoles)
        patient, _, _, _  = helpers.getPatientAndPlan(basePlanName)
        sys.stdout.write(f' \n\n ===================== [{patient.Name}] start for {newPlanName} ===================== \n')

        if 1:
            airOverride(setMaterial=True)

        # Step 1 - Copy Plan
        print (f' \n\n ===================== start for {newPlanName} ===================== \n')
        keyPlanCSTerm     = config.SUFFIX_PLAN_CS.format(config.PREFIX_CLINICAL_CONTOURS)
        keyPlanCSAutoTerm = config.SUFFIX_PLAN_CS.format(config.PREFIX_AUTOMATED_CONTOURS)
        createArcBeam     = False
        if keyPlanCSTerm in newPlanName or keyPlanCSAutoTerm in newPlanName:
            createArcBeam = True
        print (f' - [copyPlanAndOptimize()][{newPlanName}] createArcBeam: {createArcBeam}')
        copyPlanStatus = helpers.copyPlan(basePlanName, newPlanName, createArcBeam=createArcBeam, debug=debug)
        if not copyPlanStatus:
            return False, optimizeValue, timeTaken
        
        # Step 2 - Upload/update objectives
        objectiveStatus = uploadORUpdateObjectives(newPlanName, pathKNOObjectives, uploadObjectivesBool, updateObjectivesBool, forceObjectives, objectiveFType)
        if not objectiveStatus:
            return False, optimizeValue, timeTaken
        
        # Step 3 - Optimize
        optimizeStatus, optimizeValue = optimizePlan(newPlanName, optSteps, reset=optReset, pathIsoDoseXML=pathIsoDoseXML)
        if optimizeStatus:
            _ = helpers.rayStationSave()
        
        timeTaken = round(time.time() - t0, 2)
        print (f' \n\n ===================== end for {newPlanName} (in {timeTaken} s) ===================== \n\n')
        return optimizeStatus, optimizeValue, timeTaken
    
    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace()

    timeTaken = round(time.time() - t0, 2)
    return False, optimizeValue, timeTaken

# Func 0
def main(params):

    # Step 1.1 - Init patient details
    try:

        if 1:
            pathPatient = params[config.KEYNAME_PATH_PATIENT]
            forceUploadPatient = params[config.KEYNAME_FORCE_UPLOAD_PATIENT]
            forceCurrentPatient = params[config.KEYNAME_FORCE_CURRENT_PATIENT]

            # Step 1.2 - Init plan details
            keynameCancerType = params[config.KEYNAME_CANCER_TYPE]
            planNameOG        = config.KEYNAME_PLAN_OG.format(keynameCancerType)

            planNameCS        = config.KEYNAME_PLAN_CS.format(keynameCancerType, config.PREFIX_CLINICAL_CONTOURS)
            planNameDFO       = config.KEYNAME_PLAN_DFO.format(keynameCancerType, config.PREFIX_CLINICAL_CONTOURS)
            planNameDFO2      = config.KEYNAME_PLAN_DFO2.format(keynameCancerType, config.PREFIX_CLINICAL_CONTOURS)
            planNameEUD       = config.KEYNAME_PLAN_EUD.format(keynameCancerType, config.PREFIX_CLINICAL_CONTOURS)
            planNameFinal     = config.KEYNAME_PLAN_FINAL.format(keynameCancerType, config.PREFIX_CLINICAL_CONTOURS)
            planNameFinal2    = config.KEYNAME_PLAN_FINAL2.format(keynameCancerType, config.PREFIX_CLINICAL_CONTOURS)
            
            planNameCSAuto     = config.KEYNAME_PLAN_CS.format(keynameCancerType, config.PREFIX_AUTOMATED_CONTOURS)
            planNameDFOAuto    = config.KEYNAME_PLAN_DFO.format(keynameCancerType, config.PREFIX_AUTOMATED_CONTOURS)
            planNameDFO2Auto   = config.KEYNAME_PLAN_DFO2.format(keynameCancerType, config.PREFIX_AUTOMATED_CONTOURS)
            planNameEUDAuto    = config.KEYNAME_PLAN_EUD.format(keynameCancerType, config.PREFIX_AUTOMATED_CONTOURS)
            planNameFinalAuto  = config.KEYNAME_PLAN_FINAL.format(keynameCancerType, config.PREFIX_AUTOMATED_CONTOURS)
            planNameFinal2Auto = config.KEYNAME_PLAN_FINAL2.format(keynameCancerType, config.PREFIX_AUTOMATED_CONTOURS)
            
            optStepsRe        = params[config.KEYNAME_OPT_STEPS_RE]

            # Step 1.3 - Init objective details
            pathKNOObjectivesClassSolution = params[config.KEYNAME_PATH_CLASSSOL]
            pathKNOObjectives              = params[config.KEYNAME_PATH_OBJECTIVES]
                
            # Step 1.4 - Init evaluation details
            pathDVHParams = params[config.KEYNAME_PATH_DVHPARAMS]
            pathIsoDoseXML = params[config.KEYNAME_PATH_ISODOSEXML]

            # Step 1.5 - Get contour type
            contourTypeNow = params[config.KEYNAME_CONTOUR_TYPE]

            # Step 1.6 - Updates (for auto-contours)
            if contourTypeNow is not config.KEYNAME_CONTOUR_EVAL:
                pathKNOObjectivesClassSolutionAuto = helpers.updateKNOXMLForAutoContours(pathKNOObjectivesClassSolution, pathKNOObjectives, potentialRoisToRenameInAuto=config.PHOTON_POTENTIAL_ROIS_TO_RENAME_FOR_AUTO)
                if pathKNOObjectivesClassSolutionAuto is None:
                    return 0
                pathKNOObjectivesAuto = helpers.updateKNOXMLForAutoContours(pathKNOObjectives, pathKNOObjectives, potentialRoisToRenameInAuto=config.PHOTON_POTENTIAL_ROIS_TO_RENAME_FOR_AUTO)
                if pathKNOObjectivesAuto is None:
                    return 0

            # Step 1.6 - for console printing (while debugging across multiple consoles)
            try:
                sys.stdout.write(f' \n\n ===================== start for {pathPatient} [{planNameOG}]===================== \n')
            except:
                traceback.print_exc()

        ############## Step 2 - Upload data to RayStation
        try:

            tTotal = time.time()

            if contourTypeNow in [config.KEYNAME_CONTOUR_CLINICAL, config.KEYNAME_CONTOUR_AUTO]:
                uploadRTAppsStatus = uploadRTAppsDataToRStation(pathPatient, planName=planNameOG, forceUpload=forceUploadPatient, forceCurrentPatient=forceCurrentPatient)
                print ('\n - [main] uploadRTAppsStatus: ', uploadRTAppsStatus)
                if contourTypeNow == config.KEYNAME_CONTOUR_CLINICAL:    
                    print ('\n\n ------------------------- ')
                    print ('  - Computing dose on clinical contours')
                    print (' ------------------------- \n\n')
                    if uploadRTAppsStatus:
                        classSolPlanStatus, classSolPlanValue, classSolPlanTime = copyPlanAndOptimize(planNameOG, planNameCS
                                                    , pathKNOObjectivesClassSolution, uploadObjectivesBool=True, updateObjectivesBool=False, forceObjectives=False, objectiveFType=None
                                                    , optSteps=optStepsRe, optReset=True, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                        print (' - [main] classSolPlanStatus: ', classSolPlanStatus)
                        if classSolPlanStatus:
                            dfoPlanStatus, dfoPlanValue, dfoPlanTime = copyPlanAndOptimize(planNameCS, planNameDFO
                                                    , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] dfoPlanStatus: ', dfoPlanStatus)
                            if dfoPlanStatus:
                                dfo2PlanStatus, dfo2PlanValue, dfo2PlanTime = copyPlanAndOptimize(planNameDFO, planNameDFO2
                                                        , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                print (' - [main] dfo2PlanStatus: ', dfo2PlanStatus)
                                if dfo2PlanStatus:
                                    eudPlanStatus, eudPlanValue, eudPlanTime = copyPlanAndOptimize(planNameDFO2, planNameEUD
                                                        , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                    print (' - [main] eudPlanStatus: ', eudPlanStatus)
                                    if eudPlanStatus:
                                        finalPlanStatus, finalPlanValue, finalPlanTime = copyPlanAndOptimize(planNameEUD, planNameFinal
                                                            , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                        print (' - [main] finalPlanStatus: ', finalPlanStatus)
                                        if finalPlanStatus:
                                            res = helpers.evaluatePlans(pathDVHParams
                                                        , planNames=[planNameOG, planNameCS, planNameDFO, planNameDFO2, planNameEUD, planNameFinal, planNameFinalAuto]
                                                        , planTimes ={planNameOG:-1, planNameCS:classSolPlanTime, planNameDFO:dfoPlanTime, planNameDFO2:dfo2PlanTime, planNameEUD:eudPlanTime, planNameFinal:finalPlanTime, planNameFinalAuto:-1}
                                                        , planValues={planNameOG:-1, planNameCS:classSolPlanValue, planNameDFO:dfoPlanValue, planNameDFO2:dfo2PlanValue, planNameEUD:eudPlanValue, planNameFinal:finalPlanValue, planNameFinalAuto:-1}
                                                        , pathPatient=pathPatient
                                                        , contourType=config.KEYNAME_CONTOUR_CLINICAL
                                                        )

                elif contourTypeNow == config.KEYNAME_CONTOUR_AUTO:
                    print ('\n\n ------------------------- ')
                    print ('  - Computing dose on auto contours')
                    print (' ------------------------- \n\n')
                    if uploadRTAppsStatus:
                        autoContouringStatus, autoContouringTime = doAutoContouring()
                        print ('\n - [main] autoContouringStatus: ', autoContouringStatus)
                        if autoContouringStatus:
                            classSolPlanStatusAuto, classSolPlanValueAuto, classSolPlanTimeAuto = copyPlanAndOptimize(planNameOG, planNameCSAuto
                                                        , pathKNOObjectivesClassSolutionAuto, uploadObjectivesBool=True, updateObjectivesBool=False, forceObjectives=False, objectiveFType=None
                                                        , optSteps=optStepsRe, optReset=True, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] classSolPlanStatusAuto: ', classSolPlanStatusAuto)
                            if classSolPlanStatusAuto:
                                dfoPlanStatusAuto, dfoPlanValueAuto, dfoPlanTimeAuto = copyPlanAndOptimize(planNameCSAuto, planNameDFOAuto
                                                        , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                print (' - [main] dfoPlanStatusAuto: ', dfoPlanStatusAuto)
                                if dfoPlanStatusAuto:
                                    dfo2PlanStatusAuto, dfo2PlanValueAuto, dfo2PlanTimeAuto = copyPlanAndOptimize(planNameDFOAuto, planNameDFO2Auto
                                                            , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                    print (' - [main] dfo2PlanStatusAuto: ', dfo2PlanStatusAuto)
                                    if dfo2PlanStatusAuto:
                                        eudPlanStatusAuto, eudPlanValueAuto, eudPlanTimeAuto = copyPlanAndOptimize(planNameDFO2Auto, planNameEUDAuto
                                                            , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                        print (' - [main] eudPlanStatusAuto: ', eudPlanStatusAuto)
                                        if eudPlanStatusAuto:
                                            finalPlanStatusAuto, finalPlanValueAuto, finalPlanTimeAuto = copyPlanAndOptimize(planNameEUDAuto, planNameFinalAuto
                                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                            print (' - [main] finalPlanStatuAutos: ', finalPlanStatusAuto)
                                            if finalPlanStatusAuto:
                                                resAuto = helpers.evaluatePlans(pathDVHParams
                                                            , planNames=[planNameOG, planNameCSAuto, planNameDFOAuto, planNameDFO2Auto, planNameEUDAuto, planNameFinalAuto, planNameFinal]
                                                            , planTimes  = {planNameOG:-1, planNameCSAuto: classSolPlanTimeAuto , planNameDFOAuto: dfoPlanTimeAuto , planNameDFO2Auto:dfo2PlanTimeAuto , planNameEUDAuto: eudPlanTimeAuto , planNameFinalAuto: finalPlanTimeAuto , planNameFinal:-1}
                                                            , planValues = {planNameOG:-1, planNameCSAuto: classSolPlanValueAuto, planNameDFOAuto: dfoPlanValueAuto, planNameDFO2Auto:dfo2PlanValueAuto, planNameEUDAuto: eudPlanValueAuto, planNameFinalAuto: finalPlanValueAuto, planNameFinal:-1}
                                                            , planExtras = {config.KEYNAME_AUTOCONTOURING_TIME:autoContouringTime}
                                                            , pathPatient=pathPatient
                                                            , contourType=config.KEYNAME_CONTOUR_AUTO)
                                                if DEBUG_PDB: pdb.set_trace()

            elif contourTypeNow == config.KEYNAME_CONTOUR_ALL:
                
                uploadRTAppsStatus = uploadRTAppsDataToRStation(pathPatient, planName=planNameOG, forceUpload=forceUploadPatient, forceCurrentPatient=forceCurrentPatient)
                print ('\n - [main] uploadRTAppsStatus: ', uploadRTAppsStatus)

                if 1:    
                    print ('\n\n ------------------------- ')
                    print ('  - Computing dose on clinical contours')
                    print (' ------------------------- \n\n')
                    if uploadRTAppsStatus:
                        classSolPlanStatus, classSolPlanValue, classSolPlanTime = copyPlanAndOptimize(planNameOG, planNameCS
                                                    , pathKNOObjectivesClassSolution, uploadObjectivesBool=True, updateObjectivesBool=False, forceObjectives=False, objectiveFType=None
                                                    , optSteps=optStepsRe, optReset=True, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                        print (' - [main] classSolPlanStatus: ', classSolPlanStatus)
                        if classSolPlanStatus:
                            dfoPlanStatus, dfoPlanValue, dfoPlanTime = copyPlanAndOptimize(planNameCS, planNameDFO
                                                    , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] dfoPlanStatus: ', dfoPlanStatus)
                            if dfoPlanStatus:
                                dfo2PlanStatus, dfo2PlanValue, dfo2PlanTime = copyPlanAndOptimize(planNameDFO, planNameDFO2
                                                        , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                print (' - [main] dfo2PlanStatus: ', dfo2PlanStatus)
                                if dfo2PlanStatus:
                                    eudPlanStatus, eudPlanValue, eudPlanTime = copyPlanAndOptimize(planNameDFO2, planNameEUD
                                                        , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                    print (' - [main] eudPlanStatus: ', eudPlanStatus)
                                    if eudPlanStatus:
                                        finalPlanStatus, finalPlanValue, finalPlanTime = copyPlanAndOptimize(planNameEUD, planNameFinal
                                                            , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                        print (' - [main] finalPlanStatus: ', finalPlanStatus)
                                        if finalPlanStatus:
                                            res = helpers.evaluatePlans(pathDVHParams
                                                        , planNames=[planNameOG, planNameCS, planNameDFO, planNameDFO2, planNameEUD, planNameFinal]
                                                        , planTimes ={planNameOG:-1, planNameCS:classSolPlanTime, planNameDFO:dfoPlanTime, planNameDFO2:dfo2PlanTime, planNameEUD:eudPlanTime, planNameFinal:finalPlanTime}
                                                        , planValues={planNameOG:-1, planNameCS:classSolPlanValue, planNameDFO:dfoPlanValue, planNameDFO2:dfo2PlanValue, planNameEUD:eudPlanValue, planNameFinal:finalPlanValue}
                                                        , pathPatient=pathPatient
                                                        , contourType=config.KEYNAME_CONTOUR_CLINICAL)
                                            print ('\n - [main] Total time passed: ', round(time.time()-tTotal, 2), ' seconds')
                
                if 1:
                    print ('\n\n ------------------------- ')
                    print ('  - Computing dose on auto contours')
                    print (' ------------------------- \n\n')
                    if uploadRTAppsStatus:
                        autoContouringStatus, autoContouringTime = doAutoContouring()
                        print ('\n - [main] autoContouringStatus: ', autoContouringStatus)
                        if autoContouringStatus:
                            classSolPlanStatusAuto, classSolPlanValueAuto, classSolPlanTimeAuto = copyPlanAndOptimize(planNameOG, planNameCSAuto
                                                        , pathKNOObjectivesClassSolutionAuto, uploadObjectivesBool=True, updateObjectivesBool=False, forceObjectives=False, objectiveFType=None
                                                        , optSteps=optStepsRe, optReset=True, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] classSolPlanStatusAuto: ', classSolPlanStatusAuto)
                            if classSolPlanStatusAuto:
                                dfoPlanStatusAuto, dfoPlanValueAuto, dfoPlanTimeAuto = copyPlanAndOptimize(planNameCSAuto, planNameDFOAuto
                                                        , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                print (' - [main] dfoPlanStatusAuto: ', dfoPlanStatusAuto)
                                if dfoPlanStatusAuto:
                                    dfo2PlanStatusAuto, dfo2PlanValueAuto, dfo2PlanTimeAuto = copyPlanAndOptimize(planNameDFOAuto, planNameDFO2Auto
                                                            , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                    print (' - [main] dfo2PlanStatusAuto: ', dfo2PlanStatusAuto)
                                    if dfo2PlanStatusAuto:
                                        eudPlanStatusAuto, eudPlanValueAuto, eudPlanTimeAuto = copyPlanAndOptimize(planNameDFO2Auto, planNameEUDAuto
                                                            , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                        print (' - [main] eudPlanStatusAuto: ', eudPlanStatusAuto)
                                        if eudPlanStatusAuto:
                                            finalPlanStatusAuto, finalPlanValueAuto, finalPlanTimeAuto = copyPlanAndOptimize(planNameEUDAuto, planNameFinalAuto
                                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                            print (' - [main] finalPlanStatuAutos: ', finalPlanStatusAuto)
                                            if finalPlanStatusAuto:
                                                resAuto = helpers.evaluatePlans(pathDVHParams
                                                            , planNames=[planNameOG, planNameCSAuto, planNameDFOAuto, planNameDFO2Auto, planNameEUDAuto, planNameFinalAuto]
                                                            , planTimes  = {planNameOG:-1, planNameCSAuto: classSolPlanTimeAuto , planNameDFOAuto: dfoPlanTimeAuto , planNameDFO2Auto:dfo2PlanTimeAuto , planNameEUDAuto: eudPlanTimeAuto , planNameFinalAuto: finalPlanTimeAuto}
                                                            , planValues = {planNameOG:-1, planNameCSAuto: classSolPlanValueAuto, planNameDFOAuto: dfoPlanValueAuto, planNameDFO2Auto:dfo2PlanValueAuto, planNameEUDAuto: eudPlanValueAuto, planNameFinalAuto: finalPlanValueAuto}
                                                            , planExtras = {config.KEYNAME_AUTOCONTOURING_TIME:autoContouringTime}
                                                            , pathPatient=pathPatient
                                                            , contourType=config.KEYNAME_CONTOUR_AUTO)
                                                print ('\n - [main] Total time passed: ', round(time.time()-tTotal, 2), ' seconds')

                if 1:
                    if uploadRTAppsStatus:
                        print ('\n\n ------------------------- ')
                        print ('  - Compare clinical stats ')
                        print (' ------------------------- \n\n')
                        res = helpers.evaluatePlans(pathDVHParams
                                , planNames=[planNameOG, planNameCS, planNameCSAuto, planNameFinal, planNameFinalAuto]
                                , pathPatient=pathPatient
                                , contourType=config.KEYNAME_CONTOUR_ALL)

                        print ('\n\n ------------------------- ')
                        print ('  - NTCP Modeling ')
                        print (' ------------------------- \n\n')
                        patientID = Path(pathPatient).parts[-2]
                        plansForNTCP = [planNameOG, planNameFinal, planNameFinalAuto]
                        getNTCPVals(patientID, plansForNTCP, pathPatient)

                        print ('\n - [main] Total time passed: ', round(time.time()-tTotal, 2), ' seconds')
                        
            elif contourTypeNow == config.KEYNAME_CONTOUR_EVAL:
                res = helpers.evaluatePlans(pathDVHParams
                            , planNames=[planNameOG, planNameCS, planNameCSAuto, planNameFinal, planNameFinalAuto]
                            , pathPatient=pathPatient
                            , contourType=config.KEYNAME_CONTOUR_ALL)

                print ('\n\n ------------------------- ')
                print ('  - NTCP Modeling ')
                print (' ------------------------- \n\n')
                patientID = Path(pathPatient).parts[-2]
                plansForNTCP = [planNameOG, planNameFinal, planNameFinalAuto]
                getNTCPVals(patientID, plansForNTCP, pathPatient)

                print ('\n - [main] Total time passed: ', round(time.time()-tTotal, 2), ' seconds')

            elif contourTypeNow == config.KEYNAME_CONTOUR_DEBUG:
                
                # upload
                if 0:
                    uploadRTAppsStatus = uploadRTAppsDataToRStation(pathPatient, planName=planNameOG, forceUpload=forceUploadPatient, forceCurrentPatient=forceCurrentPatient)
                    print ('\n - [main] uploadRTAppsStatus: ', uploadRTAppsStatus)

                # for R2-R5 updates
                if 0:
                    try:
                        dfoPlanStatus, dfoPlanValue, dfoPlanTime = copyPlanAndOptimize(planNameCS, planNameDFO
                                                    , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                        print (' - [main] dfoPlanStatus: ', dfoPlanStatus)
                        if dfoPlanStatus:
                            dfo2PlanStatus, dfo2PlanValue, dfo2PlanTime = copyPlanAndOptimize(planNameDFO, planNameDFO2
                                                    , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] dfo2PlanStatus: ', dfo2PlanStatus)
                            if dfo2PlanStatus:
                                eudPlanStatus, eudPlanValue, eudPlanTime = copyPlanAndOptimize(planNameDFO2, planNameEUD
                                                    , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                print (' - [main] eudPlanStatus: ', eudPlanStatus)
                                if eudPlanStatus:
                                    finalPlanStatus, finalPlanValue, finalPlanTime = copyPlanAndOptimize(planNameEUD, planNameFinal
                                                        , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                    print (' - [main] finalPlanStatus: ', finalPlanStatus)
                                    if finalPlanStatus:
                                        res = helpers.evaluatePlans(pathDVHParams
                                                    , planNames=[planNameOG, planNameCS, planNameDFO, planNameDFO2, planNameEUD, planNameFinal, planNameFinalAuto]
                                                    , planTimes ={planNameOG:-1, planNameDFO:dfoPlanTime, planNameDFO2:dfo2PlanTime, planNameEUD:eudPlanTime, planNameFinal:finalPlanTime}
                                                    , planValues={planNameOG:-1, planNameDFO:dfoPlanValue, planNameDFO2:dfo2PlanValue, planNameEUD:eudPlanValue, planNameFinal:finalPlanValue}
                                                    , pathPatient=pathPatient
                                                    , contourType=config.KEYNAME_CONTOUR_CLINICAL
                                                    , save=True
                                                    )
                                        print ('\n - [main] Total time passed: ', round(time.time()-tTotal, 2), ' seconds')
                    except:
                        traceback.print_exc()

                # upload + auto-contouring
                if 0:
                    uploadRTAppsStatus = uploadRTAppsDataToRStation(pathPatient, planName=planNameOG, forceUpload=forceUploadPatient, forceCurrentPatient=forceCurrentPatient)
                    print ('\n - [main] uploadRTAppsStatus: ', uploadRTAppsStatus)
                    if uploadRTAppsStatus:
                        autoContouringStatus, autoContouringTime = doAutoContouring()
                        print ('\n - [main] autoContouringStatus: ', autoContouringStatus)

                # for A1-A5 updates
                if 0:
                    autoContouringStatus, autoContouringTime = doAutoContouring()
                    print ('\n - [main] autoContouringStatus: ', autoContouringStatus)
                    if autoContouringStatus:
                        classSolPlanStatusAuto, classSolPlanValueAuto, classSolPlanTimeAuto = copyPlanAndOptimize(planNameOG, planNameCSAuto
                                                    , pathKNOObjectivesClassSolutionAuto, uploadObjectivesBool=True, updateObjectivesBool=False, forceObjectives=False, objectiveFType=None
                                                    , optSteps=optStepsRe, optReset=True, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                        print (' - [main] classSolPlanStatusAuto: ', classSolPlanStatusAuto)
                        if classSolPlanStatusAuto:
                            dfoPlanStatusAuto, dfoPlanValueAuto, dfoPlanTimeAuto = copyPlanAndOptimize(planNameCSAuto, planNameDFOAuto
                                                    , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] dfoPlanStatusAuto: ', dfoPlanStatusAuto)
                            if dfoPlanStatusAuto:
                                dfo2PlanStatusAuto, dfo2PlanValueAuto, dfo2PlanTimeAuto = copyPlanAndOptimize(planNameDFOAuto, planNameDFO2Auto
                                                        , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                print (' - [main] dfo2PlanStatusAuto: ', dfo2PlanStatusAuto)
                                if dfo2PlanStatusAuto:
                                    eudPlanStatusAuto, eudPlanValueAuto, eudPlanTimeAuto = copyPlanAndOptimize(planNameDFO2Auto, planNameEUDAuto
                                                        , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                    print (' - [main] eudPlanStatusAuto: ', eudPlanStatusAuto)
                                    if eudPlanStatusAuto:
                                        finalPlanStatusAuto, finalPlanValueAuto, finalPlanTimeAuto = copyPlanAndOptimize(planNameEUDAuto, planNameFinalAuto
                                                            , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                        print (' - [main] finalPlanStatuAutos: ', finalPlanStatusAuto)
                                        if finalPlanStatusAuto:
                                            resAuto = helpers.evaluatePlans(pathDVHParams
                                                        , planNames=[planNameOG, planNameCSAuto, planNameDFOAuto, planNameDFO2Auto, planNameEUDAuto, planNameFinalAuto, planNameFinal]
                                                        , planTimes  = {planNameOG:-1, planNameCSAuto: classSolPlanTimeAuto , planNameDFOAuto: dfoPlanTimeAuto , planNameDFO2Auto:dfo2PlanTimeAuto , planNameEUDAuto: eudPlanTimeAuto , planNameFinalAuto: finalPlanTimeAuto , planNameFinal:-1}
                                                        , planValues = {planNameOG:-1, planNameCSAuto: classSolPlanValueAuto, planNameDFOAuto: dfoPlanValueAuto, planNameDFO2Auto:dfo2PlanValueAuto, planNameEUDAuto: eudPlanValueAuto, planNameFinalAuto: finalPlanValueAuto, planNameFinal:-1}
                                                        , planExtras = {config.KEYNAME_AUTOCONTOURING_TIME:autoContouringTime}
                                                        , pathPatient=pathPatient
                                                        , contourType=config.KEYNAME_CONTOUR_AUTO)
                                            if DEBUG_PDB: pdb.set_trace()
                
                # For A2/A3/A4/A5 updates
                if 0:
                    try:
                        dfoPlanStatusAuto, dfoPlanValueAuto, dfoPlanTimeAuto = copyPlanAndOptimize(planNameCSAuto, planNameDFOAuto
                                                    , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                        print (' - [main] dfoPlanStatusAuto: ', dfoPlanStatusAuto)
                        if dfoPlanStatusAuto:
                            dfo2PlanStatusAuto, dfo2PlanValueAuto, dfo2PlanTimeAuto = copyPlanAndOptimize(planNameDFOAuto, planNameDFO2Auto
                                                    , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] dfo2PlanStatusAuto: ', dfo2PlanStatusAuto)
                            if dfo2PlanStatusAuto:
                                eudPlanStatusAuto, eudPlanValueAuto, eudPlanTimeAuto = copyPlanAndOptimize(planNameDFO2Auto, planNameEUDAuto
                                                    , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                print (' - [main] eudPlanStatusAuto: ', eudPlanStatusAuto)
                                if eudPlanStatusAuto:
                                    finalPlanStatusAuto, finalPlanValueAuto, finalPlanTimeAuto = copyPlanAndOptimize(planNameEUDAuto, planNameFinalAuto
                                                        , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                    print (' - [main] finalPlanStatuAutos: ', finalPlanStatusAuto)
                                    if finalPlanStatusAuto:
                                        print ('\n\n ------------------------- ')
                                        print ('  - Compare clinical stats ')
                                        print (' ------------------------- \n\n')
                                        res = helpers.evaluatePlans(pathDVHParams
                                            , planNames=[planNameOG, planNameCS, planNameCSAuto, planNameFinal, planNameFinalAuto]
                                            , pathPatient=pathPatient
                                            , contourType=config.KEYNAME_CONTOUR_ALL
                                        )

                                        print ('\n\n ------------------------- ')
                                        print ('  - NTCP Modeling ')
                                        print (' ------------------------- \n\n')
                                        patientID = Path(pathPatient).parts[-2]
                                        plansForNTCP = [planNameOG, planNameFinal, planNameFinalAuto]
                                        getNTCPVals(patientID, plansForNTCP, pathPatient)

                                        print ('\n - [main] Total time passed: ', round(time.time()-tTotal, 2), ' seconds')
                    except:
                        traceback.print_exc()
                
                # For A3/A4/A5 updates
                if 0:
                    try:
                        dfo2PlanStatusAuto, dfo2PlanValueAuto, dfo2PlanTimeAuto = copyPlanAndOptimize(planNameDFOAuto, planNameDFO2Auto
                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_DOSEFALLOFF
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                        print (' - [main] dfo2PlanStatusAuto: ', dfo2PlanStatusAuto)
                        if dfo2PlanStatusAuto:
                            eudPlanStatusAuto, eudPlanValueAuto, eudPlanTimeAuto = copyPlanAndOptimize(planNameDFO2Auto, planNameEUDAuto
                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] eudPlanStatusAuto: ', eudPlanStatusAuto)
                            if eudPlanStatusAuto:
                                finalPlanStatusAuto, finalPlanValueAuto, finalPlanTimeAuto = copyPlanAndOptimize(planNameEUDAuto, planNameFinalAuto
                                                    , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                print (' - [main] finalPlanStatuAutos: ', finalPlanStatusAuto)
                                if finalPlanStatusAuto:
                                    print ('\n\n ------------------------- ')
                                    print ('  - Compare clinical stats ')
                                    print (' ------------------------- \n\n')
                                    res = helpers.evaluatePlans(pathDVHParams
                                        , planNames=[planNameOG, planNameCS, planNameCSAuto, planNameFinal, planNameFinalAuto]
                                        , pathPatient=pathPatient
                                        , contourType=config.KEYNAME_CONTOUR_ALL
                                    )

                                    print ('\n\n ------------------------- ')
                                    print ('  - NTCP Modeling ')
                                    print (' ------------------------- \n\n')
                                    patientID = Path(pathPatient).parts[-2]
                                    plansForNTCP = [planNameOG, planNameFinal, planNameFinalAuto]
                                    getNTCPVals(patientID, plansForNTCP, pathPatient)

                                    print ('\n - [main] Total time passed: ', round(time.time()-tTotal, 2), ' seconds')
                    except:
                        traceback.print_exc()

                # for A4/A5 updates
                if 0:
                    try:
                        eudPlanStatusAuto, eudPlanValueAuto, eudPlanTimeAuto = copyPlanAndOptimize(planNameDFO2Auto, planNameEUDAuto
                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                        print (' - [main] eudPlanStatusAuto: ', eudPlanStatusAuto)
                        if eudPlanStatusAuto:
                            finalPlanStatusAuto, finalPlanValueAuto, finalPlanTimeAuto = copyPlanAndOptimize(planNameEUDAuto, planNameFinalAuto
                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] finalPlanStatusAuto: ', finalPlanStatusAuto)
                            if finalPlanStatusAuto:
                                print ('\n\n ------------------------- ')
                                print ('  - Compare clinical stats ')
                                print (' ------------------------- \n\n')
                                res = helpers.evaluatePlans(pathDVHParams
                                    , planNames=[planNameOG, planNameCS, planNameCSAuto, planNameFinal, planNameFinalAuto]
                                    , pathPatient=pathPatient
                                    , contourType=config.KEYNAME_CONTOUR_ALL
                                )

                                print ('\n\n ------------------------- ')
                                print ('  - NTCP Modeling ')
                                print (' ------------------------- \n\n')
                                patientID = Path(pathPatient).parts[-2]
                                plansForNTCP = [planNameOG, planNameFinal, planNameFinalAuto]
                                getNTCPVals(patientID, plansForNTCP, pathPatient)

                                print ('\n - [main] Total time passed: ', round(time.time()-tTotal, 2), ' seconds')

                    except:
                        traceback.print_exc()

                # NTCP
                if 0:
                    patientID = Path(pathPatient).parts[-2]
                    patientRSObj = helpers.loadPatientUsingID(patientID)
                    if patientRSObj is not None:
                        plansForNTCP = [planNameOG, planNameFinal, planNameFinalAuto]
                        getNTCPVals(patientID, plansForNTCP, pathPatient)
                    else:
                        print (' - [main] patientRSObj is None')

        except:
            traceback.print_exc()
            if DEBUG_PDB: pdb.set_trace()
    
    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace()

if __name__ == "__main__":

    ###################################################################################
    # Step 0 - Init project paths
    ###################################################################################
    DIR_THIS = Path(__file__).parent.absolute() # P:\RayStationScripts
    DIR_DATA = Path(DIR_THIS).parent.absolute().joinpath('RayStationData')
    
    ###################################################################################
    # Step 1 - Logging details
    ###################################################################################
    if 1:
        DIR_LOGS        = Path(DIR_THIS).joinpath('_logs', 'logsPhoton', 'run6')
        loggerTimestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        pathLogFile     = Path(DIR_LOGS).joinpath("photon-log_{}.txt".format(loggerTimestamp))
        Path(pathLogFile).parent.mkdir(parents=True, exist_ok=True)
        logging.basicConfig(filename=str(pathLogFile), level=logging.DEBUG, filemode='a', format='%(asctime)s[%(levelname)s] %(funcName)s: %(message)s',datefmt='%d/%m/%Y %I:%M:%S %p')

    ###################################################################################
    # Step 2 - Get common files
    ###################################################################################
    pathKNOObjectivesClassSolution = Path(DIR_DATA).joinpath('assets', 'objective-template-photon-kno.xml')
    pathDVHParams = Path(DIR_DATA).joinpath('assets', 'eval-template-photon.csv')
    pathIsoDoseXML = Path(DIR_DATA).joinpath('assets', 'isodose.xml'); # pathIsoDoseXML = None

    ###################################################################################
    # Step 3 -  Specific patient paths (of data extracted from RTPACS)
    ###################################################################################

    if 1:
        # HCAI-Dose-1
        if 1:
            pathPatient = Path(DIR_DATA).joinpath('')
            pathKNOObjectivesClinical = Path(pathPatient).joinpath('.xml')
            
            keynameCancerType = '1A OROFARKL'
            optStepsRe        = 4

            forceUploadPatient   = True
            forceCurrentPatient  = False

            contourType = config.KEYNAME_CONTOUR_ALL # config.KEYNAME_CONTOUR_EVAL, config.KEYNAME_CONTOUR_AUTO, config.KEYNAME_CONTOUR_CLINICAL, config.KEYNAME_CONTOUR_ALL
            pathPatient = Path(DIR_DATA).joinpath('LUMC-Dose', 'HCAI-Dose-x40', '2.25.143009416517888739377675690506069860146')
            pathKNOObjectivesClinical = Path(pathPatient).joinpath('1017729.xml')

            keynameCancerType = '1A HYPOFARKL'
            optStepsRe        = 4

            forceUploadPatient   = False
            forceCurrentPatient  = True

            contourType = config.KEYNAME_CONTOUR_EVAL # [config.{KEYNAME_CONTOUR_EVAL, KEYNAME_CONTOUR_AUTO, KEYNAME_CONTOUR_CLINICAL, KEYNAME_CONTOUR_ALL, KEYNAME_CONTOUR_DEBUG}]

        params = {

            # Patient Case related
            config.KEYNAME_PATH_PATIENT           : pathPatient
            , config.KEYNAME_FORCE_UPLOAD_PATIENT : forceUploadPatient
            , config.KEYNAME_FORCE_CURRENT_PATIENT: forceCurrentPatient
        
            # Paths
            , config.KEYNAME_PATH_CLASSSOL        : pathKNOObjectivesClassSolution
            , config.KEYNAME_PATH_OBJECTIVES      : pathKNOObjectivesClinical
            , config.KEYNAME_PATH_DVHPARAMS       : pathDVHParams
            , config.KEYNAME_PATH_ISODOSEXML      : pathIsoDoseXML
            
            # Plan parameters
            , config.KEYNAME_CANCER_TYPE            : keynameCancerType
            , config.KEYNAME_OPT_STEPS_RE           : optStepsRe
            , config.KEYNAME_CONTOUR_TYPE           : contourType
        }
        print ('\n -------------------------- [main] params: ')
        print (params)

        ###################################################################################
        # Step 3 - Run main
        ###################################################################################
        main(params)

# To run code in RS console (and print to console)
"""
import sys
import runpy
from pathlib import Path

sys.path.append(str(Path('H:\\').joinpath('RayStationScripts')))
pathFile = str(Path('H:\\').joinpath('RayStationScripts', 'hnDoseEval-v4.py'))
_ = runpy.run_path(pathFile, run_name='__main__')

"""