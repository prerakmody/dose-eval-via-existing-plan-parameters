"""
This script is used for proton planning (using beam settings and objective function weights from an existing proton plan)
It expects
    - the clinical plan to already be present in RStation (along with objectives)
    - assets/objective-template-proton-kno.xml
    - assets/eval-template-proton-robust.json
    - assets/isodose.xml

Run the script from "__main__"
"""

# Import Raystation libraries
import connect

# Import private libraries
import hnDoseEvalHelpers as helpers
import hnDoseConfig as config

# Import public libraries
import re
import sys
import pdb
import time
import json
import shutil
import logging
import datetime
import traceback
import numpy as np
from pathlib import Path

DEBUG_PDB = True
def print(*args, **kwargs):
    logging.info(" ".join(map(str, args)), **kwargs)

########################################################
#                        HELPERS                       #
########################################################

def checkIfRoiIgnore(planName, objectiveFType, roiName):
    
    # Step 1 - Init
    roisFullIgnore = []
    roisPartialIgnore = []
    roiIgnoreBool = False

    planEUDSuffix      = config.SUFFIX_PLAN_DFO.format(config.PREFIX_CLINICAL_CONTOURS) # R2
    planEUDAutoSuffix  = config.SUFFIX_PLAN_DFO.format(config.PREFIX_AUTOMATED_CONTOURS)
    planEUD2Suffix       = config.SUFFIX_PLAN_DFO2.format(config.PREFIX_CLINICAL_CONTOURS) # R3
    planEUD2AutoSuffix   = config.SUFFIX_PLAN_DFO2.format(config.PREFIX_AUTOMATED_CONTOURS)
    planFinalSuffix     = config.SUFFIX_PLAN_FINAL.format(config.PREFIX_CLINICAL_CONTOURS)  # R5
    planFinalAutoSuffix = config.SUFFIX_PLAN_FINAL.format(config.PREFIX_AUTOMATED_CONTOURS)

    # Step 2 - Ignore check on certain ROIs
    if roiName in [config.KEYNAME_BRAINSTEM, config.KEYNAME_BRAINSTEM + config.KEY_3MM_SUFFIX
                   , config.KEYNAME_SPINALCORD, config.KEYNAME_SPINALCORD + config.KEY_3MM_SUFFIX
                   , config.KEYNAME_GHOST_CRANIAL, config.KEYNAME_EAR_L_GHOST, config.KEYNAME_EAR_R_GHOST
                   ]:
        return roiIgnoreBool

    # Step 2 - Get ROIs to ignore (while updating objectives)
    roisPartialIgnoreBasic = [
                config.REGEX_ROI_HOT, config.REGEX_ROI_HEET
                , config.REGEX_KOUD, config.REGEX_COLD
                , config.REGEX_ROI_MINDOSE, config.REGEX_ROI_MAXDOSE
                , config.REGEX_ONDER, config.REGEX_OVER
                , config.REGEX_DFO
                , config.REGEX_ONDERDOSERING
                , config.REGEX_ROI_HOSTPSOTS
                , config.REGEX_CAUDAAL, config.REGEX_CAUD
                ]
    if (planEUDSuffix in planName or planEUDAutoSuffix in planName) and objectiveFType == config.KEY_FTYPE_MAXEUD:
        roisFullIgnore = [config.KEYNAME_BODY] # does not matter as it is Opt_Body in protons
        roisPartialIgnore = roisPartialIgnoreBasic + [config.KEYNAME_3MM] # the extra part compared to EUD2 is KEYNAME_3MM
    
    elif (planEUD2Suffix in planName or planEUD2AutoSuffix in planName) and objectiveFType == config.KEY_FTYPE_MAXEUD:
        roisFullIgnore = [config.KEYNAME_BODY, config.KEYNAME_RING_PTV_DL2, config.KEYNAME_GHOST]
        roisPartialIgnore = roisPartialIgnoreBasic
    
    elif (planFinalSuffix in planName or planFinalAutoSuffix in planName) and objectiveFType == None:
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

def updateObjectivesForProton(planName, objectivesFromRS, objectivesFromPath, objectiveFType):
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
        planDFO2TermSuffix      = config.SUFFIX_PLAN_DFO2.format(config.PREFIX_CLINICAL_CONTOURS)
        planDFO2AutoTermSuffix  = config.SUFFIX_PLAN_DFO2.format(config.PREFIX_AUTOMATED_CONTOURS)
        planFinalTermSuffix     = config.SUFFIX_PLAN_FINAL.format(config.PREFIX_CLINICAL_CONTOURS)  # R5
        planFinalAutoTermSuffix = config.SUFFIX_PLAN_FINAL.format(config.PREFIX_AUTOMATED_CONTOURS) # A5

        # Step 2 - Get objectives from file (a.k.a path)
        updatedFromPathObjectivesStatus = {}
        for pathIdx, objectiveFromPath in enumerate(objectivesFromPath):
            newRoiName      = objectiveFromPath.roi_name
            newFunctionType = objectiveFromPath.function_type
            newWeight       = objectiveFromPath.weight
            newRobust       = objectiveFromPath.is_robust
            if newWeight == 0:
                continue
            
            newKey = f'{newRoiName}-{newFunctionType}-{newRobust}'
            updatedFromPathObjectivesStatus[newKey] = {'exists': False, 'idx':pathIdx}

        # Step 4 - Get existing objectives and update them
        if len(objectivesFromRS) > 0:
            print ('\n\n - [INFO][updateObjectivesForProton()][Patient={}] Updating Objectives ... '.format(helpers.getPatientIdentifier(patient)))
            
            # Step 4.1 - Update weights of existing objectives (based on function-type and roi-name)
            for rsIdx, objectiveFromRS in enumerate(objectivesFromRS):
                try:

                    # Step 4.1.1 - Get roiName, weight, fType
                    existingRoiName      = objectiveFromRS.ForRegionOfInterest.Name
                    existingWeight       = objectiveFromRS.DoseFunctionParameters.Weight
                    existingRobust       = objectiveFromRS.UseRobustness 
                    try   : existingFunctionType = objectiveFromRS.DoseFunctionParameters.FunctionType
                    except: existingFunctionType = config.KEY_FTYPE_DOSEFALLOFF
                    if existingRoiName == config.KEYNAME_BODY: continue
                    # [DEBUG]
                    if 0:
                        print (' - [DEBUG][updateObjectivesForProton()] \troi: {}, \tfunctionType: {}, \tweight: {}'.format(existingRoiName, existingFunctionType, existingWeight))
                    
                    # Step 4.1.2 - Make note of existing objectives that are in new objectives
                    existingKey = f'{existingRoiName}-{existingFunctionType}-{existingRobust}'
                    if existingKey in updatedFromPathObjectivesStatus:
                        updatedFromPathObjectivesStatus[existingKey]['exists'] = True

                    # Step 4.1.3 - Check if existing objective is in new objectives
                    for objectiveFromPath in objectivesFromPath:
                        newRoiName      = objectiveFromPath.roi_name
                        newFunctionType = objectiveFromPath.function_type
                        newWeight       = objectiveFromPath.weight
                        newRobust       = objectiveFromPath.is_robust
                        if newWeight == 0: 
                            continue
                        functionMatchBool = newFunctionType == objectiveFType
                        if  planFinalTermSuffix in planName or planFinalAutoTermSuffix in planName:
                            functionMatchBool = True # i.e update all rois
                        
                        robustMatchBool = existingRobust == newRobust

                        if functionMatchBool and robustMatchBool:
                            if not checkIfRoiIgnore(planName, objectiveFType, newRoiName):
                                if existingRoiName == newRoiName and existingFunctionType == newFunctionType:

                                    # Step 4.1.1 - Update weight 
                                    objectivesFromRS[rsIdx].DoseFunctionParameters.Weight = newWeight
                                    print (f'  -- [updateObjectivesForProton()] roi: {existingRoiName}, \tfType: {existingFunctionType}, \tweight: {existingWeight} --> {newWeight}')

                                    # Step 4.1.2 - Update dose level(s)
                                    if existingFunctionType == config.KEY_FTYPE_MAXEUD:
                                        existingDoseLevel = objectivesFromRS[rsIdx].DoseFunctionParameters.DoseLevel
                                        existingEudParameterA = objectivesFromRS[rsIdx].DoseFunctionParameters.EudParameterA
                                        objectivesFromRS[rsIdx].DoseFunctionParameters.DoseLevel     = objectiveFromPath.doselevel
                                        objectivesFromRS[rsIdx].DoseFunctionParameters.EudParameterA = objectiveFromPath.eud_parameter_a
                                        print (f' --- [updateObjectivesForProton()] \tDoseLevel: {existingDoseLevel} --> {objectiveFromPath.doselevel}')
                                        print (f' --- [updateObjectivesForProton()] \tEudParameterA: {existingEudParameterA} --> {objectiveFromPath.eud_parameter_a}')
                                    elif existingFunctionType == config.KEY_FTYPE_DOSEFALLOFF:
                                        if planDFO2TermSuffix in planName or planDFO2AutoTermSuffix in planName: # [Ref: https://iprova.lumc.nl/Portal/#/document/1fc74366-40c1-4b84-a653-0455b6c891f8 (Vervolgens voor beide opties)]
                                            print (f' --- [updateObjectivesForProton()] \tHighDoseLevel: {objectiveFromRS.DoseFunctionParameters.HighDoseLevel} --> {objectiveFromPath.high_doselevel}')
                                            print (f' --- [updateObjectivesForProton()] \tLowDoseLevel : {objectiveFromRS.DoseFunctionParameters.LowDoseLevel}  --> {objectiveFromPath.low_doselevel}')
                                            objectivesFromRS[rsIdx].DoseFunctionParameters.HighDoseLevel = objectiveFromPath.high_doselevel
                                            objectivesFromRS[rsIdx].DoseFunctionParameters.LowDoseLevel  = objectiveFromPath.low_doselevel
                                    
                                    elif existingFunctionType == config.KEY_FTYPE_MAXDOSE:
                                        pass
                                
                                # Step 4.1.3 - For debugging purposes
                                key = f'{newRoiName}-{newFunctionType}-{newRobust}'
                                updatedFromPathObjectivesStatus[key]['updated'] = True
                except:
                    traceback.print_exc()
                    if DEBUG_PDB: pdb.set_trace()

        else:
            print (f' - [updateKNOObjectivesinRStation()] No objectives found in RS for plan {planName}')

        
        # Step 5 - Add new objectives [if 1) they dont exist, and if they match 2) function-type and 3) roiName conditions]
        print ('')
        for key in updatedFromPathObjectivesStatus:
            if updatedFromPathObjectivesStatus[key]['exists'] is False:
                roiName = key.split('-')[0]
                fType   = key.split('-')[1]
                isRobust = key.split('-')[2]
                
                fTypeCondition     = fType == objectiveFType # only if ftType in [DoseFallOff, MaxEUD]
                roiIgnoreCondition = checkIfRoiIgnore(planName, objectiveFType, roiName) # this roiName comes from existing objectives of RS
                if objectiveFType == None: # since this is the last plan, we want to add all (non-existent) objectives
                    fTypeCondition = True

                if fTypeCondition and not roiIgnoreCondition:
                    print (f' - [updateObjectivesForProton()] Adding new objective for {roiName} of type {fType}')
                    newObjective = objectivesFromPath[updatedFromPathObjectivesStatus[key]['idx']]
                    print (f' --- [updateObjectivesForProton()] newObjective for roi: {newObjective.roi_name} and fType: {newObjective.function_type} with weight: {newObjective.weight}')
                    if newObjective.weight > 0:
                        newObjective.apply()

        updateObjectiveStatus = True

    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace()
    

    return updateObjectiveStatus

def uploadOrUpdateProtonObjectives(planName, pathKNOProtonObjectives, uploadObjectivesBool, updateObjectivesBool, forceObjectives, objectiveFType):

    # Step 1 - Init
    objectiveStatus    = False

    # Step 2 - Boilerplate code for get patient and plans
    _, case, plan, beamset  = helpers.getPatientAndPlan(planName)

    # Step 3 - Get objectives
    objectivesFromPath = helpers.getObjectivesFromPath(plan, beamset, pathKNOProtonObjectives)
    objectivesFromRS   = helpers.getObjectivesFromPlan(plan, beamset)

    # Step 3.1 - Filter an objective
    if 1:
        for objIdx, objectiveFromPath in enumerate(objectivesFromPath):
            try:
                if objectiveFromPath.roi_name == config.KEYNAME_OPT_BODY and objectiveFromPath.function_type == config.KEY_FTYPE_MAXDOSE:
                    if float(objectiveFromPath.doselevel) < 7000:
                        del objectivesFromPath[objIdx]
            except:
                pass

    # Step 4 - Upload/Update Objectives
    if uploadObjectivesBool:
        if len(objectivesFromRS) == 0 or forceObjectives:
            objectiveStatus = helpers.uploadObjectivesToRS(plan, beamset, objectivesFromPath)
        else:
            print(' - [uploadORUpdateObjectives()] Objectives already exist, not uploading')
    
    elif updateObjectivesBool:
        if len(objectivesFromRS) > 0:
            objectiveStatus = updateObjectivesForProton(planName, objectivesFromRS, objectivesFromPath, objectiveFType)
        else:
            print(' - [uploadORUpdateObjectives()] No objectives found, not updating')
    
    # Step 5 - Do sanity check
    helpers.doPlanSanityCheck(case, planName)

    return objectiveStatus

def robustEvaluationOld(planName, force=False):

    try:
        
        # Step 0 - Init
        patientObj, caseObj, planObj, beamsetObj  = helpers.getPatientAndPlan(planName)

        if force:
            radiationSetScenarioGroups = caseObj.TreatmentDelivery.RadiationSetScenarioGroups
            for id_ in range(len(radiationSetScenarioGroups)):
                if radiationSetScenarioGroups[id_].Name == planName:
                    print (f' - [robustEvaluation()] Deleting scenario group: {planName}')
                    radiationSetScenarioGroups[id_].DeleteRadiationSetScenarioGroup()
                    helpers.rayStationSave()
        
        boolScenarioAvailable = False
        for id_ in range(len(radiationSetScenarioGroups)):
            if radiationSetScenarioGroups[id_].Name == planName:
                boolScenarioAvailable = True
                break

        if boolScenarioAvailable is False:
            print (f' - [robustEvaluation()] Creating scenario group: {planName}')
            retval_0 = beamsetObj.CreateRadiationSetScenarioGroup(Name=planName
                    , UseIsotropicPositionUncertainty=False
                    , PositionUncertaintySuperior=0.3, PositionUncertaintyInferior=0.3
                    , PositionUncertaintyPosterior=0.3, PositionUncertaintyAnterior=0.3
                    , PositionUncertaintyLeft=0.3, PositionUncertaintyRight=0.3
                    , PositionUncertaintyFormation="AxesAndDiagonalEndPoints", PositionUncertaintyList=None
                    , DensityUncertainty=3, NumberOfDensityDiscretizationPoints=2
                    , ComputeScenarioDosesAfterGroupCreation=False)

            print (f' - [robustEvaluation()] Computing dose using IonMonteCarlo algo')
            beamsetObj.ComputeDose(ComputeBeamDoses=True, DoseAlgorithm="IonMonteCarlo", ForceRecompute=True)

            print (f' - [robustEvaluation()] Computing scenario group dose values: {planName}')
            retval_0.ComputeScenarioGroupDoseValues()
            helpers.rayStationSave()
            print (f' - [robustEvaluation()] Done computing scenario group dose values: {planName}')
        
        else:
            print (f' - [robustEvaluation()] Scenario group already exists: {planName}')

    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace()

def robustEvaluationViaUI(planName, force=False):

    try:
        
        # Step 1 - Init
        _, caseObj, _ , beamSetObj = helpers.getPatientAndPlan(planName, debug=True)

        # Step 2 - Delete scenario group if force
        if force:
            radiationSetScenarioGroups = caseObj.TreatmentDelivery.RadiationSetScenarioGroups
            for id_ in range(len(radiationSetScenarioGroups)):
                try:
                    if radiationSetScenarioGroups[id_].Name == planName:
                        print (f' - [robustEvaluationViaUI()] Deleting scenario group: {planName}')
                        radiationSetScenarioGroups[id_].DeleteRadiationSetScenarioGroup()
                        helpers.rayStationSave()
                        break
                except:
                    pass

        # Step 3 - Confirm scenarioGroup does not exist
        boolScenarioAvailable = False
        radiationSetScenarioGroups = caseObj.TreatmentDelivery.RadiationSetScenarioGroups
        for id_ in range(len(radiationSetScenarioGroups)):
            if radiationSetScenarioGroups[id_].Name == planName:
                boolScenarioAvailable = True
                break
        
        # Step 4 - Create scenario group if not exists (and do robust eval)
        if boolScenarioAvailable is False:
            print (f' - [robustEvaluationViaUI()] Creating scenario group: {planName}')
            retval_0 = beamSetObj.CreateRadiationSetScenarioGroup(Name=planName
                    , UseIsotropicPositionUncertainty=False
                    , PositionUncertaintySuperior=0.3, PositionUncertaintyInferior=0.3
                    , PositionUncertaintyPosterior=0.3, PositionUncertaintyAnterior=0.3
                    , PositionUncertaintyLeft=0.3, PositionUncertaintyRight=0.3
                    , PositionUncertaintyFormation="AxesAndDiagonalEndPoints", PositionUncertaintyList=None
                    , DensityUncertainty=3, NumberOfDensityDiscretizationPoints=2
                    , ComputeScenarioDosesAfterGroupCreation=False)

            print (f' - [robustEvaluationViaUI()] Computing dose using IonMonteCarlo algo')
            beamSetObj.ComputeDose(ComputeBeamDoses=True, DoseAlgorithm="IonMonteCarlo", ForceRecompute=True)

            print (f' - [robustEvaluationViaUI()] Computing scenario group dose values: {planName}')
            retval_0.ComputeScenarioGroupDoseValues()
            helpers.rayStationSave()
            print (f' - [robustEvaluationViaUI()] Done computing scenario group dose values: {planName}')
        else:
            print (f' - [robustEvaluationViaUI()] Scenario group already exists: {planName}')

    except:
        traceback.print_exc()

def robustEvaluationViaSelf(planName, pathRobustTemplate, pathRobustResultsSave, verbose=False):

    try:
        
        def getClinicalGoalsFromDose(doseEvalObj):
            """
            For evalFunc in planObj.TreatmentCourse.EvaluationSetup.EvaluationFunctions:
                - evalFunc.ForRegionOfInterest.Name    : self-explanatory
                - evalFunc.PlanningGoal.Type           : 'DoseAtVolume', 'DoseAtAbsoluteVolume', 'DoseAtRelativeVolume', 'VolumeAtDose', 'AbsoluteVolumeAtDose', 'RelativeVolumeAtDose'
                - evalFunc.PlanningGoal.GoalCriteria   : 'AtMost', 'AtLeast'
                - evalFunc.PlanningGoal.AcceptanceLevel: float
            """
            
            # Step 0 - Init
            resObj = {}

            try:
                
                # Step 1 - Init
                _, caseObj, planObj, beamSetObj = helpers.getPatientAndPlan(planName)
                evalFuncs = planObj.TreatmentCourse.EvaluationSetup.EvaluationFunctions

                # Step 2 - Loop over funcs
                for evalFunc in evalFuncs:
                    try:
                        value = evalFunc.GetClinicalGoalValueForEvaluationDose(DoseDistribution=doseEvalObj, ScaleFractionDoseToBeamSet=True)
                    except:
                        value = evalFunc.GetClinicalGoalValue()
                    roiGoalkey   = '-'.join([evalFunc.ForRegionOfInterest.Name, evalFunc.PlanningGoal.Type, evalFunc.PlanningGoal.GoalCriteria, str(evalFunc.PlanningGoal.AcceptanceLevel), str(evalFunc.PlanningGoal.ParameterValue)])
                    tmp          = [value, evalFunc.PlanningGoal.GoalCriteria, evalFunc.PlanningGoal.AcceptanceLevel]
                    if evalFunc.PlanningGoal.AcceptanceLevel > 0:
                        if evalFunc.PlanningGoal.GoalCriteria == config.KEYNAME_ATMOST:  
                            if value <= evalFunc.PlanningGoal.AcceptanceLevel: resObj[roiGoalkey] = tmp + [config.KEYNAME_PASS]
                            else                                             : resObj[roiGoalkey] = tmp + [config.KEYNAME_FAIL]
                        
                        elif evalFunc.PlanningGoal.GoalCriteria == config.KEYNAME_ATLEAST:
                            if value >= evalFunc.PlanningGoal.AcceptanceLevel: resObj[roiGoalkey] = tmp + [config.KEYNAME_PASS]
                            else                                             : resObj[roiGoalkey] = tmp + [config.KEYNAME_FAIL]

            except:
                traceback.print_exc()
                if DEBUG_PDB: pdb.set_trace()
            
            return resObj


        def getKeyMapping(robustTemplateObj, thisObj):
            
            # Step 0 - Init
            mapping = {}

            try:
                
                # Step 1 - Loop over keys
                roiGoalValueKeys = list(thisObj.keys())
                for roiGoalValueKey in roiGoalValueKeys:
                    thisRoiName         = roiGoalValueKey.split('-')[0]
                    thisGoalType        = roiGoalValueKey.split('-')[1]
                    thisGoalCriteria    = roiGoalValueKey.split('-')[2]
                    thisAcceptanceValue = roiGoalValueKey.split('-')[3]
                    thisParameterValue  = roiGoalValueKey.split('-')[4]
                    mapping[roiGoalValueKey] = ''

                    for templateKey in robustTemplateObj.keys():
                        templateKeyRoiName = templateKey.split(' ')[0]
                        templateKeyGoal    = templateKey[templateKey.find('[')+1: templateKey.find(']')]

                        roiCheckBool = thisRoiName == templateKeyRoiName
                        if roiCheckBool:
                            if thisRoiName == 'CTV_DL1':
                                perc = int(round(float(thisParameterValue)/config.CTV_DL1_MAX,2) * 100)
                                if str(perc) in templateKeyGoal:
                                    mapping[roiGoalValueKey] = templateKey
                                    break
                            elif thisRoiName == 'CTV_DL2':
                                perc = int(round(float(thisParameterValue)/config.CTV_DL2_MAX,2) * 100)
                                if str(perc) in templateKeyGoal:
                                    mapping[roiGoalValueKey] = templateKey
                                    break
                            elif thisRoiName == 'Brainstem_Core':
                                if str(int(float(thisAcceptanceValue))) in templateKeyGoal:
                                    mapping[roiGoalValueKey] = templateKey
                                    break
                            elif thisRoiName == 'Brainstem_Surf':
                                if str(int(float(thisAcceptanceValue))) in templateKeyGoal:
                                    mapping[roiGoalValueKey] = templateKey
                                    break
                            elif thisRoiName == 'SpinalCord_Core':
                                if str(int(float(thisAcceptanceValue))) in templateKeyGoal:
                                    mapping[roiGoalValueKey] = templateKey
                                    break
                            elif thisRoiName == 'SpinalCord_Surf':
                                if str(int(float(thisAcceptanceValue))) in templateKeyGoal:
                                    mapping[roiGoalValueKey] = templateKey
                                    break
                            else:
                                mapping[roiGoalValueKey] = templateKey
                                break

            except:
                traceback.print_exc()
            
            return mapping

        # Step 1 - Init
        _, caseObj, planObj, beamSetObj = helpers.getPatientAndPlan(planName, debug=True)
        radiationSetScenarioGroups = caseObj.TreatmentDelivery.RadiationSetScenarioGroups
        res       = {config.KEYNAME_NOMINAL_DOSE: {}, config.KEYNAME_MAX_DOSE: {}, config.KEYNAME_MIN_DOSE: {}, config.KEYNAME_SCENARIO_DOSE: {}}
        dose3DMin = []
        dose3DMax = []

        # Step 2 - FInd idx of radiationSetScenarioGroups
        groupIdx = -1
        for id_ in range(len(radiationSetScenarioGroups)):
            try:
                if radiationSetScenarioGroups[id_].Name == planName:
                    groupIdx = id_
            except:
                pass
        
        if groupIdx > -1:
            print (f' - [robustEvaluationViaSelf()] Found scenario group: {planName} and doing robust eval on all scenarios')
            for scenarioId, scenarioObj in enumerate(caseObj.TreatmentDelivery.RadiationSetScenarioGroups[groupIdx].DiscreteFractionDoseScenarios):
                if verbose:
                    print (f' - [robustEvaluationViaSelf()] Computing dose for scenario: {scenarioId}')
                    # print (scenarioId, scenarioObj.PerturbedDoseProperties.IsoCenterShift, scenarioObj.PerturbedDoseProperties.RelativeDensityShift)
                
                # Step 2.2 - Compute dose for scenario
                beamSetObj.ComputePerturbedDose(DensityPerturbation=scenarioObj.PerturbedDoseProperties.RelativeDensityShift
                        , PatientShift=scenarioObj.PerturbedDoseProperties.IsoCenterShift
                        , OnlyOneDosePerImageSet=True, ExaminationNames=[caseObj.Examinations[0].Name], FractionNumbers=[0]
                )

                # Step 2.3 - Get doseObj and dose values
                doseEvalObj = caseObj.TreatmentDelivery.FractionEvaluations[0].DoseOnExaminations[0].DoseEvaluations[0]
                dose3Darray = doseEvalObj.DoseValues.DoseData

                # Step 2.4 - Compute min and max dose
                if len(dose3DMin) == 0: dose3DMin = dose3Darray
                else                  : dose3DMin = np.minimum(dose3DMin, dose3Darray)
                if len(dose3DMax) == 0: dose3DMax = dose3Darray
                else                  : dose3DMax = np.maximum(dose3DMax, dose3Darray)

                # Step 2.5 - Compute clinical goals for each perturbed dose
                res[config.KEYNAME_SCENARIO_DOSE][scenarioId] = getClinicalGoalsFromDose(doseEvalObj)

                # Step 2.99 - Debug
                # if scenarioId > -1:
                #     break

            # Step 3 - Set dose values to min and max dose
            doseEvalObj.SetDoseValues(Array=dose3DMax.flatten(), CalculationInfo='Voxelwise max', DoseAlgorithm='Undefined')
            res[config.KEYNAME_MAX_DOSE] = getClinicalGoalsFromDose(doseEvalObj)

            doseEvalObj.SetDoseValues(Array=dose3DMin.flatten(), CalculationInfo='Voxelwise min', DoseAlgorithm='Undefined')
            res[config.KEYNAME_MIN_DOSE] = getClinicalGoalsFromDose(doseEvalObj)

            # Step 4 - Get norminal dose
            doseEvalObj = planObj.TreatmentCourse.TotalDose # type=CompositeDose
            res[config.KEYNAME_NOMINAL_DOSE] = getClinicalGoalsFromDose(doseEvalObj)

            # Step 5 - Compute pass/fail (and save voxelwise-man/max values)
            if Path(pathRobustTemplate).exists():
                with open(str(pathRobustTemplate), 'r') as fp:
                    robustTemplate = json.load(fp)

                    # Step 5.0 - Init
                    keyMapper = getKeyMapping(robustTemplate, res[config.KEYNAME_SCENARIO_DOSE][0])

                    # Step 5.1 - Compute pass percentage
                    roiKeys = list(res[config.KEYNAME_SCENARIO_DOSE][0].keys())
                    for roiKey in roiKeys:
                        
                        try:
                            roiTemplateKey = keyMapper[roiKey]

                            if len(roiTemplateKey) > 0:
                                
                                # Step 5.1.1 - Collect for KEYNAME_SCENARIOS
                                roiPassCount = 0
                                for scenarioId in res[config.KEYNAME_SCENARIO_DOSE].keys():
                                    if res[config.KEYNAME_SCENARIO_DOSE][scenarioId][roiKey][-1] == config.KEYNAME_PASS:
                                        roiPassCount += 1
                                    robustTemplate[roiTemplateKey][config.KEYNAME_SCENARIOS].append(round(res[config.KEYNAME_SCENARIO_DOSE][scenarioId][roiKey][0],4))
                                robustTemplate[roiTemplateKey][config.KEYNAME_PASSED] = round(roiPassCount/len(res[config.KEYNAME_SCENARIO_DOSE].keys()), 4)

                                # Step 5.1.2 - Collect for KEYNAME_NOMINAL_DOSE
                                robustTemplate[roiTemplateKey][config.KEYNAME_NOMINAL_DOSE] = round(res[config.KEYNAME_NOMINAL_DOSE][roiKey][0], 4)

                                # Step 5.1.3 - Collect for KEYNAME_MIN_DOSE and KEYNAME_MAX_DOSE
                                if roiTemplateKey in config.OBJECTIVES_ROIS_TUMOR:
                                    robustTemplate[roiTemplateKey][config.KEYNAME_VOXELWISE_WORST] = round(res[config.KEYNAME_MIN_DOSE][roiKey][0], 4)
                                else:
                                    robustTemplate[roiTemplateKey][config.KEYNAME_VOXELWISE_WORST] = round(res[config.KEYNAME_MAX_DOSE][roiKey][0], 4)
                        except:
                            traceback.print_exc()
                            if DEBUG_PDB: pdb.set_trace()
        
                    # Step 6 - Save
                    sys.stdout.write(f' - [robustEvaluationViaSelf()] Saving robustTemplate to: {pathRobustResultsSave}')
                    print (f' - [robustEvaluationViaSelf()] Saving robustTemplate to: {pathRobustResultsSave}')
                    with open(str(pathRobustResultsSave), 'w') as fp:
                        robustTemplateWrite = {planName: robustTemplate}
                        json.dump(robustTemplateWrite, fp, indent=4)
            
            else:
                print (f' - [robustEvaluationViaSelf()] Robust template not found: {pathRobustTemplate}')
        
        else:
            print (f' - [robustEvaluationViaSelf()] Scenario group not found for {planName}')

    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace()

def robustEvaluation(planName, pathRobustTemplate, pathPatient, force=False, verbose=False):

    try:
        
        t0 = time.time()
        # Step 0 - Init
        pathRobustResultsSave = Path(pathPatient) / config.FILENAME_ROBUST_EVAL_RESULTS.format(planName)

        # Step 1 - Do UI-based robust eval
        robustEvaluationViaUI(planName, force=force)
        
        # Step 2 - Do self-based robust eval
        robustEvaluationViaSelf(planName, pathRobustTemplate, pathRobustResultsSave, verbose=verbose)
        print (f' - [robustEvaluation()] Done robust eval for {planName} in {round(time.time() - t0, 2)} s')

    except:
        traceback.print_exc()
        pdb.set_trace()

########################################################
#                          NTCP                        #
########################################################

def getNTCPVals(patientID, plans, pathPatient):

    try:
        
        print (f' - [getNTCPVals()] Getting NTCP values for patient: {patientID}')
        
        # Step 0 - Plans and tumor types
        patientObjOld = {
            'HCAI-Dose-P1': ['1F PVFOTONEN', None] # None
            , 'HCAI-Dose-P2': ['1A HYPOFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P3': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P4': ['1A KNOKL', None] # None
            , 'HCAI-Dose-P5': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P6': ['1A HYPOFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P7': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P8': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P9': ['1A HYPOFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P10': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P11': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P12': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P13': ['1A LARYNXKL', helpers.NtcpKnoDysphagia.TumorLocationType.LARYNX]
            , 'HCAI-Dose-P14': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P15': ['1A LARYNXKL', helpers.NtcpKnoDysphagia.TumorLocationType.LARYNX]
            , 'HCAI-Dose-P16': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P17': ['1A LARYNXKL', helpers.NtcpKnoDysphagia.TumorLocationType.LARYNX]
            , 'HCAI-Dose-P18': ['1A KNOKL', None] # None
            , 'HCAI-Dose-P19': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P20': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P21': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P22': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P23': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P24': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P25': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P26': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P27': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P28': ['1A TONGKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P29': ['1A LARYNXKL', helpers.NtcpKnoDysphagia.TumorLocationType.LARYNX]
            , 'HCAI-Dose-P30': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P31': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
            , 'HCAI-Dose-P32': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.ORAL_CAVITY]
        }

        patientObj = {
            'HCAI-Dose-P1': ['1F PVFOTONEN', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX] # None
            , 'HCAI-Dose-P2': ['1A HYPOFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P3': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P4': ['1A KNOKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX] # None
            , 'HCAI-Dose-P5': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P6': ['1A HYPOFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P7': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P8': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P9': ['1A HYPOFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P10': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P11': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P12': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P13': ['1A LARYNXKL', helpers.NtcpKnoDysphagia.TumorLocationType.LARYNX]
            , 'HCAI-Dose-P14': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P15': ['1A LARYNXKL', helpers.NtcpKnoDysphagia.TumorLocationType.LARYNX]
            , 'HCAI-Dose-P16': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P17': ['1A LARYNXKL', helpers.NtcpKnoDysphagia.TumorLocationType.LARYNX]
            , 'HCAI-Dose-P18': ['1A KNOKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX] # None
            , 'HCAI-Dose-P19': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P20': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P21': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P22': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P23': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P24': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P25': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P26': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P27': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P28': ['1A TONGKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P29': ['1A LARYNXKL', helpers.NtcpKnoDysphagia.TumorLocationType.LARYNX]
            , 'HCAI-Dose-P30': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P31': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
            , 'HCAI-Dose-P32': ['1A OROFARKL', helpers.NtcpKnoDysphagia.TumorLocationType.PHARYNX]
        }
        
        # Step 1 - Other params
        params = {
            config.KEY_NTCP_PLAN1  : None
            , config.KEY_NTCP_PLAN2: None
            , config.KEY_NTCP_TREATMENT_TYPE     : helpers.TreatmentType.PRIMARY
            , config.KEY_NTCP_TUMOR_LOCATION     : patientObj[patientID][1] # NtcpKnoDysphagia.TumorLocationType.{ORAL_CAVITY, PHARYNX, LARYNX}
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
                objNTCP = helpers.KNONTCP(params)

                # Step 2.2 - Calculate
                tmp = objNTCP.ntcp_xerostomia_grade_2_values
                if config.KEY_NTCP_PLAN_1 in tmp:
                    res[plan][tmp[config.KEY_NTCP_MODEL_NAME]] = round(tmp[config.KEY_NTCP_PLAN_1][config.KEY_NTCP_PLAN_VALUE], 5)

                tmp = objNTCP.ntcp_xerostomia_grade_3_values
                if config.KEY_NTCP_PLAN_1 in tmp:
                    res[plan][tmp[config.KEY_NTCP_MODEL_NAME]] = round(tmp[config.KEY_NTCP_PLAN_1][config.KEY_NTCP_PLAN_VALUE], 5)
                
                tmp = objNTCP.ntcp_dysphagia_grade_2_values
                if config.KEY_NTCP_PLAN_1 in tmp:
                    res[plan][tmp[config.KEY_NTCP_MODEL_NAME]] = round(tmp[config.KEY_NTCP_PLAN_1][config.KEY_NTCP_PLAN_VALUE], 5)
                
                tmp = objNTCP.ntcp_dysphagia_grade_3_values
                if config.KEY_NTCP_PLAN_1 in tmp:
                    res[plan][tmp[config.KEY_NTCP_MODEL_NAME]] = round(tmp[config.KEY_NTCP_PLAN_1][config.KEY_NTCP_PLAN_VALUE], 5)
                
        # Step 3 - Save
        pathPatientNTCP = Path(pathPatient) / config.FILENAME_NTCP_RESULTS
        with open(str(pathPatientNTCP), 'w') as fp:
            json.dump(res, fp, indent=4)
        
        print (f' - [getNTCPVals()] Done NTCP calculations for {patientID}')
        
    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace()

########################################################
#                    AUTO-HELPERS                      #
########################################################

def doAutoContouringForProton():

    autoContourStatus = False
    t0 = time.time()
    try:
        
        # Step 0 - Init
        print (f' \n\n ===================== start autocontouring (for protons) ===================== \n\n')
        patient = helpers.rayStationSave()
        case    = patient.Cases[0]
        examinationName = case.Examinations[0].Name

        # Step 1.1 - Get OARs to include    
        oarDuplicateStatus = helpers.checkOARDuplicateStatus(case)
        oarsToInclude = []
        for oar in oarDuplicateStatus:
            if not oarDuplicateStatus[oar]:
                oarsToInclude.append(oar)
        print (' - [doAutoContouringForProton()] OARs pending to auto-contour: ', oarsToInclude)
                
        # Step 1.2 - Run OAR segmentation
        if len(oarsToInclude):
            _ = helpers.rayStationSave()
            case.SetCurrent() # Why == Potential Error: "The case to which the examination belongs mist be selected"
            examination = case.Examinations[0]
            _ = examination.RunOarSegmentation(ModelName="RSL Head and Neck CT", ExaminationsAndRegistrations={ examinationName: None }, RoisToInclude=oarsToInclude)
            helpers.rayStationSave()
        else:
            print (' - [doAutoContouringForProton()] No OARs to auto-contour')
                
        # Step 1.3 - Check aut--contours status
        ## If they dont have the config.KEY_AUTOCONTOUR_SUFFIX, then they were not initially present, so simply rename them
        oarDuplicateStatus = helpers.checkOARDuplicateStatus(case)
        for oar in oarDuplicateStatus:
            if not oarDuplicateStatus[oar]:
                print (' - [doAutoContouringForProton()] OAR auto-contour not present: ', oar)
                if oar == config.KEYNAME_CAVITY_ORAL:
                    newName                                            = config.KEYNAME_ORAL_CAVITY + config.KEY_AUTOCONTOUR_SUFFIX
                    case.PatientModel.RegionsOfInterest[oar].Name      = newName
                    case.PatientModel.RegionsOfInterest[newName].Color = config.OAR_DUPLICATE_COLOR_RGB_STRING 
                    print (' - [doAutoContouringForProton()] Renamed OAR: ', oar, ' --> ', newName)
                elif oar == config.KEYNAME_ESOPHAGUS_S:
                    try:
                        newName                                            = config.KEYNAME_ESOPHAGUS + config.KEY_AUTOCONTOUR_SUFFIX
                        case.PatientModel.RegionsOfInterest[oar].Name      = newName
                        case.PatientModel.RegionsOfInterest[newName].Color = config.OAR_DUPLICATE_COLOR_RGB_STRING 
                        print (' - [doAutoContouringForProton()] Renamed OAR: ', oar, ' --> ', newName)
                    except:
                        traceback.print_exc()
                else:
                    case.PatientModel.RegionsOfInterest[oar].Name = oar + config.KEY_AUTOCONTOUR_SUFFIX

        # Step 1.4 - Do ROI Algebra
        helpers.doROIAlgebraForProtonAutoContours(case)
        _ = helpers.rayStationSave()
        timeTaken = round(time.time() - t0, 2)
        print (f' \n\n ===================== end autocontouring (for protons) (in {timeTaken} s) ===================== \n\n')
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

        # if len(pathPatientCTFolders) == 1 and len(pathPatientRTDoseFolders) == 1 and len(pathPatientRTPlanFolders) == 1 and len(pathPatientRTStructFolders) == 1:
        if len(pathPatientCTFolders) == 1 and len(pathPatientRTStructFolders) == 1:
            pathPatientCTFolder       = pathPatientCTFolders[0]
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
                        if len(pathPatientRTDoseFolders) == 1 and len(pathPatientRTPlanFolders) == 1:
                            pathPatientRTDoseFolder   = pathPatientRTDoseFolders[0]
                            pathPatientRTPlanFolder   = pathPatientRTPlanFolders[0]
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
                            print ('\n\n [uploadRTAppsDataToRStation()][Patient={}] --------------------- Step 1.3: No RTPlan/RTDose .dcm found \n\n'.format(helpers.getPatientIdentifier(patient)))
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
    
    if len(pathPatientRTDoseFolders) == 1 and len(pathPatientRTPlanFolders) == 1:
        if patientCTBool and patientRTStructBool and patientRTPlanBool and patientRTDoseBool:
            return True
        else:
            return False
    else:
        if patientCTBool and patientRTStructBool:
            return True
        else:
            return False

# Func 2
def copyProtonPlanAndOptimize(basePlanName, newPlanName
                        , pathKNOProtonObjectives, uploadObjectivesBool, updateObjectivesBool, forceObjectives, objectiveFType
                        , optSteps, optReset, pathIsoDoseXML=None
                        , debug=False):

    t0 = time.time()
    optimizeValue = -1
    timeTaken = -1
    try:

        # Step 0 - for console printing (while debugging across multiple consoles)
        try:
            patient, _, _, _  = helpers.getPatientAndPlan(basePlanName)
            sys.stdout.write(f' \n\n ===================== [{patient.Name}] start for {newPlanName} ===================== \n')
        except:
            traceback.print_exc()

        # Step 1 - Copy Plan
        print (f' \n\n ===================== start for {newPlanName} ===================== \n')
        copyPlanStatus = helpers.copyPlan(basePlanName, newPlanName, createArcBeam=False, debug=debug)
        if not copyPlanStatus:
            return False, optimizeValue, timeTaken
        
        # Step 2 - Upload/update objectives
        objectiveStatus = uploadOrUpdateProtonObjectives(newPlanName, pathKNOProtonObjectives, uploadObjectivesBool, updateObjectivesBool, forceObjectives, objectiveFType)
        if not objectiveStatus:
            return False, optimizeValue, timeTaken
        
        # Step 3 - Optimize
        if optSteps:
            optimizeStatus, optimizeValue = helpers.optimizePlan(newPlanName, optSteps, reset=optReset, pathIsoDoseXML=pathIsoDoseXML)
            if optimizeStatus:
                _ = helpers.rayStationSave()
        else:
            optimizeStatus, optimizeValue = True, -1

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

    try:
        tPatient = time.time()
        # Step 1.1 - Init patient details
        pathPatient       = params[config.KEYNAME_PATH_PATIENT]
        patientID         = params[config.KEY_PATIENTID]
        loadPatient      = params[config.KEYNAME_FORCE_LOAD_PATIENT]
        uploadPatient    = params[config.KEYNAME_FORCE_UPLOAD_PATIENT]

        # Step 1.2 - Init plan details
        keynameCancerType = params[config.KEYNAME_CANCER_TYPE]
        planNameOG        = config.KEYNAME_PLAN_OG.format(keynameCancerType)

        planNameCS        = config.KEYNAME_PLAN_CS.format(keynameCancerType, config.PREFIX_CLINICAL_CONTOURS)  # -R1
        planNameEUD       = config.KEYNAME_PLAN_DFO.format(keynameCancerType, config.PREFIX_CLINICAL_CONTOURS) # -R2
        planNameEUD2      = config.KEYNAME_PLAN_DFO2.format(keynameCancerType, config.PREFIX_CLINICAL_CONTOURS) # -R3
        planNameFinal     = config.KEYNAME_PLAN_FINAL.format(keynameCancerType, config.PREFIX_CLINICAL_CONTOURS) # -R5
        
        planNameCSAuto     = config.KEYNAME_PLAN_CS.format(keynameCancerType, config.PREFIX_AUTOMATED_CONTOURS)
        planNameEUDAuto    = config.KEYNAME_PLAN_DFO.format(keynameCancerType, config.PREFIX_AUTOMATED_CONTOURS)
        planNameEUD2Auto   = config.KEYNAME_PLAN_DFO2.format(keynameCancerType, config.PREFIX_AUTOMATED_CONTOURS)
        planNameFinalAuto  = config.KEYNAME_PLAN_FINAL.format(keynameCancerType, config.PREFIX_AUTOMATED_CONTOURS)
        
        optStepsRe        = params[config.KEYNAME_OPT_STEPS_RE]

        # Step 1.3 - Init objective details
        pathKNOObjectivesClassSolution     = params[config.KEYNAME_PATH_CLASSSOL]
        pathKNOObjectives                  = params[config.KEYNAME_PATH_OBJECTIVES]
        
        # Step 1.4 - Init evaluation details
        pathRobustTemplate = params[config.KEYNAME_PATH_ROBUST_TEMPLATE]
        pathIsoDoseXML     = params[config.KEYNAME_PATH_ISODOSEXML]
        if Path(pathIsoDoseXML).exists() is False:
            sys.stdout.write(f' - [main()] pathIsoDoseXML does not exist: {pathIsoDoseXML}')
            print (f' - [main()] pathIsoDoseXML does not exist: {pathIsoDoseXML}')
            return 0
        if Path(pathRobustTemplate).exists() is False:
            sys.stdout.write(f' - [main()] pathRobustTemplate does not exist: {pathRobustTemplate}')
            print (f' - [main()] pathRobustTemplate does not exist: {pathRobustTemplate}')
            return 0

        # Step 1.5 - Get contour type
        contourTypeNow = params[config.KEYNAME_CONTOUR_TYPE]

        # Step 1.6 - Updates (for auto-contours)
        if contourTypeNow != config.KEYNAME_CONTOUR_EVAL:
            pathKNOObjectivesClassSolutionAuto = helpers.updateKNOXMLForAutoContours(pathKNOObjectivesClassSolution, pathKNOObjectives, potentialRoisToRenameInAuto=config.PROTON_POTENTIAL_ROIS_TO_RENAME_FOR_AUTO)
            if pathKNOObjectivesClassSolutionAuto is None:
                return 0
            pathKNOObjectivesAuto = helpers.updateKNOXMLForAutoContours(pathKNOObjectives, pathKNOObjectives, potentialRoisToRenameInAuto=config.PROTON_POTENTIAL_ROIS_TO_RENAME_FOR_AUTO)
            if pathKNOObjectivesAuto is None:
                return 0

        # Step 1.7 - for console printing (while debugging across multiple consoles)
        sys.stdout.write(f' \n\n ===================== start for {pathPatient} [{planNameOG}]===================== \n')
        
        ############## Step 2 - Upload data to RayStation ##############
        try:
            
            # Step 2.1 - Get patient (for protons patients are already in RStation. Simply load them)
            if loadPatient:
                patientRSObj = helpers.loadPatientUsingID(patientID)
            elif uploadPatient:
                uploadRTAppsStatus = uploadRTAppsDataToRStation(pathPatient, planName=planNameOG, forceUpload=True, forceCurrentPatient=False)
                if uploadRTAppsStatus:
                    patientRSObj = connect.get_current(config.KEYNAME_PATIENT)
                else:
                    patientRSObj = None
            else:
                patientRSObj = connect.get_current(config.KEYNAME_PATIENT)

            if patientRSObj is None:
                sys.stdout.write(f' - [main()] Patient {patientID} not found in RayStation')
                print (f' - [main()] Patient {patientID} not found in RayStation')
                return 0

            ###################################
            if contourTypeNow in [config.KEYNAME_CONTOUR_CLINICAL, config.KEYNAME_CONTOUR_AUTO]:
                if contourTypeNow == config.KEYNAME_CONTOUR_CLINICAL:    
                    print ('\n\n ------------------------- ')
                    print ('  - Computing dose on clinical contours')
                    print (' ------------------------- \n\n')
                    classSolPlanStatus, classSolPlanValue, classSolPlanTime = copyProtonPlanAndOptimize(planNameOG, planNameCS
                                                , pathKNOObjectivesClassSolution, uploadObjectivesBool=True, updateObjectivesBool=False, forceObjectives=True, objectiveFType=None
                                                , optSteps=optStepsRe, optReset=True, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                    print (' - [main] classSolPlanStatus: ', classSolPlanStatus)
                    if classSolPlanStatus:
                        eudPlanStatus, eudPlanValue, eudPlanTime = copyProtonPlanAndOptimize(planNameCS, planNameEUD
                                        , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                        print (' - [main] eudPlanStatus: ', eudPlanStatus)
                        if eudPlanStatus:
                            eud2PlanStatus, eud2PlanValue, eud2PlanTime = copyProtonPlanAndOptimize(planNameEUD, planNameEUD2
                                                , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                            print (' - [main] eud2PlanStatus: ', eud2PlanStatus)
                            if eud2PlanStatus:
                                finalPlanStatus, finalPlanValue, finalPlanTime = copyProtonPlanAndOptimize(planNameEUD2, planNameFinal
                                                    , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                                print (' - [main] finalPlanStatus: ', finalPlanStatus)
                                if finalPlanStatus:
                                    robustEvaluation(planNameFinal, pathRobustTemplate, pathPatient, force=True)

                elif contourTypeNow == config.KEYNAME_CONTOUR_AUTO:
                    print ('\n\n ------------------------- ')
                    print ('  - Computing dose on auto contours')
                    print (' ------------------------- \n\n')
                    autoContouringStatus, autoContouringTime = doAutoContouringForProton()
                    print ('\n - [main] autoContouringStatus: ', autoContouringStatus)
                    if autoContouringStatus:
                        classSolPlanStatus, classSolPlanValue, classSolPlanTime = copyProtonPlanAndOptimize(planNameOG, planNameCSAuto
                                                    , pathKNOObjectivesClassSolutionAuto, uploadObjectivesBool=True, updateObjectivesBool=False, forceObjectives=True, objectiveFType=None
                                                    , optSteps=optStepsRe, optReset=True, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                        print (' - [main] classSolPlanStatus: ', classSolPlanStatus)
                        if classSolPlanStatus:
                            eudPlanStatus, eudPlanValue, eudPlanTime = copyProtonPlanAndOptimize(planNameCSAuto, planNameEUDAuto
                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] eudPlanStatus: ', eudPlanStatus)
                            if eudPlanStatus:
                                eud2PlanStatus, eud2PlanValue, eud2PlanTime = copyProtonPlanAndOptimize(planNameEUDAuto, planNameEUD2Auto
                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                                print (' - [main] eud2PlanStatus: ', eud2PlanStatus)
                                if eud2PlanStatus:
                                    finalPlanStatus, finalPlanValue, finalPlanTime = copyProtonPlanAndOptimize(planNameEUD2Auto, planNameFinalAuto
                                                        , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                    print (' - [main] finalPlanStatus: ', finalPlanStatus)
                                    if finalPlanStatus:
                                        robustEvaluation(planNameFinalAuto, pathRobustTemplate, pathPatient, force=True)
            
            elif contourTypeNow == config.KEYNAME_CONTOUR_ALL:
                if 1:    
                    print ('\n\n ------------------------- ')
                    print ('  - Computing dose on clinical contours')
                    print (' ------------------------- \n\n')
                    classSolPlanStatus, classSolPlanValue, classSolPlanTime = copyProtonPlanAndOptimize(planNameOG, planNameCS
                                                , pathKNOObjectivesClassSolution, uploadObjectivesBool=True, updateObjectivesBool=False, forceObjectives=True, objectiveFType=None
                                                , optSteps=optStepsRe, optReset=True, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                    print (' - [main] classSolPlanStatus: ', classSolPlanStatus)
                    if classSolPlanStatus:
                        eudPlanStatus, eudPlanValue, eudPlanTime = copyProtonPlanAndOptimize(planNameCS, planNameEUD
                                        , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                        print (' - [main] eudPlanStatus: ', eudPlanStatus)
                        if eudPlanStatus:
                            eud2PlanStatus, eud2PlanValue, eud2PlanTime = copyProtonPlanAndOptimize(planNameEUD, planNameEUD2
                                                , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                            print (' - [main] eud2PlanStatus: ', eud2PlanStatus)
                            if eud2PlanStatus:
                                finalPlanStatus, finalPlanValue, finalPlanTime = copyProtonPlanAndOptimize(planNameEUD2, planNameFinal
                                                    , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                    , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                                print (' - [main] finalPlanStatus: ', finalPlanStatus)
                                if finalPlanStatus:
                                    robustEvaluation(planNameFinal, pathRobustTemplate, pathPatient, force=True)
                                    robustEvaluation(planNameOG, pathRobustTemplate, pathPatient, force=False, verbose=True)
                
                if 1:
                    print ('\n\n ------------------------- ')
                    print ('  - Computing dose on auto contours')
                    print (' ------------------------- \n\n')
                    autoContouringStatus, autoContouringTime = doAutoContouringForProton()
                    print ('\n - [main] autoContouringStatus: ', autoContouringStatus)
                    if autoContouringStatus:
                        classSolPlanStatus, classSolPlanValue, classSolPlanTime = copyProtonPlanAndOptimize(planNameOG, planNameCSAuto
                                                    , pathKNOObjectivesClassSolutionAuto, uploadObjectivesBool=True, updateObjectivesBool=False, forceObjectives=True, objectiveFType=None
                                                    , optSteps=optStepsRe, optReset=True, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                        print (' - [main] classSolPlanStatus: ', classSolPlanStatus)
                        if classSolPlanStatus:
                            eudPlanStatus, eudPlanValue, eudPlanTime = copyProtonPlanAndOptimize(planNameCSAuto, planNameEUDAuto
                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                            print (' - [main] eudPlanStatus: ', eudPlanStatus)
                            if eudPlanStatus:
                                eud2PlanStatus, eud2PlanValue, eud2PlanTime = copyProtonPlanAndOptimize(planNameEUDAuto, planNameEUD2Auto
                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                                print (' - [main] eud2PlanStatus: ', eud2PlanStatus)
                                if eud2PlanStatus:
                                    finalPlanStatus, finalPlanValue, finalPlanTime = copyProtonPlanAndOptimize(planNameEUD2Auto, planNameFinalAuto
                                                        , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                                    print (' - [main] finalPlanStatus: ', finalPlanStatus)
                                    if finalPlanStatus:
                                        robustEvaluation(planNameFinalAuto, pathRobustTemplate, pathPatient, force=True)

                if 1:
                    plansForNTCP = [planNameOG, planNameFinal, planNameFinalAuto]
                    getNTCPVals(patientID, plansForNTCP, pathPatient)

            elif contourTypeNow == config.KEYNAME_CONTOUR_EVAL:
                robustEvaluation(planNameFinal, pathRobustTemplate, pathPatient, force=True, verbose=True)
                robustEvaluation(planNameFinalAuto, pathRobustTemplate, pathPatient, force=True, verbose=True)
                robustEvaluation(planNameOG, pathRobustTemplate, pathPatient, force=False, verbose=True)

            elif contourTypeNow == config.KEYNAME_CONTOUR_DEBUG:
                
                # upload + auto-contouring
                if 1:
                    autoContouringStatus, autoContouringTime = doAutoContouringForProton()
                    print ('\n - [main] autoContouringStatus: ', autoContouringStatus)

                # From -R2 onwards
                if 0:
                    eudPlanStatus, eudPlanValue, eudPlanTime = copyProtonPlanAndOptimize(planNameCS, planNameEUD
                                            , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                    print (' - [main] eudPlanStatus: ', eudPlanStatus)
                    if eudPlanStatus:
                        eud2PlanStatus, eud2PlanValue, eud2PlanTime = copyProtonPlanAndOptimize(planNameEUD, planNameEUD2
                                            , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                        print (' - [main] eud2PlanStatus: ', eud2PlanStatus)
                        if eud2PlanStatus:
                            finalPlanStatus, finalPlanValue, finalPlanTime = copyProtonPlanAndOptimize(planNameEUD2, planNameFinal
                                                , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                            print (' - [main] finalPlanStatus: ', finalPlanStatus)
                            if finalPlanStatus:
                                robustEvaluation(planNameFinal, pathRobustTemplate, pathPatient, force=True, verbose=True)

                # From -R3 onwards
                if 0:
                    eud2PlanStatus, eud2PlanValue, eud2PlanTime = copyProtonPlanAndOptimize(planNameEUD, planNameEUD2
                                            , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                    print (' - [main] eud2PlanStatus: ', eud2PlanStatus)
                    if eud2PlanStatus:
                        finalPlanStatus, finalPlanValue, finalPlanTime = copyProtonPlanAndOptimize(planNameEUD2, planNameFinal
                                            , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                        print (' - [main] finalPlanStatus: ', finalPlanStatus)
                        if finalPlanStatus:
                            robustEvaluation(planNameFinal, pathRobustTemplate, pathPatient, force=True, verbose=True)
                
                # From -R5 onwards
                if 0:
                    finalPlanStatus, finalPlanValue, finalPlanTime = copyProtonPlanAndOptimize(planNameEUD2, planNameFinal
                                        , pathKNOObjectives, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                    print (' - [main] finalPlanStatus: ', finalPlanStatus)
                    if finalPlanStatus:
                        robustEvaluation(planNameFinal, pathRobustTemplate, pathPatient, force=True, verbose=True)
                
                # From -R3 onwards
                if 0:
                    eud2PlanStatus, eud2PlanValue, eud2PlanTime = copyProtonPlanAndOptimize(planNameEUDAuto, planNameEUD2Auto
                                                , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=config.KEY_FTYPE_MAXEUD
                                                , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=True)
                    print (' - [main] eud2PlanStatus: ', eud2PlanStatus)
                    if eud2PlanStatus:
                        finalPlanStatus, finalPlanValue, finalPlanTime = copyProtonPlanAndOptimize(planNameEUD2Auto, planNameFinalAuto
                                            , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                            , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                        print (' - [main] finalPlanStatus: ', finalPlanStatus)
                        if finalPlanStatus:
                            robustEvaluation(planNameFinalAuto, pathRobustTemplate, pathPatient, force=True, verbose=True)

                # From -A5 onwards
                if 0:
                    finalPlanStatus, finalPlanValue, finalPlanTime = copyProtonPlanAndOptimize(planNameEUD2Auto, planNameFinalAuto
                                        , pathKNOObjectivesAuto, uploadObjectivesBool=False, updateObjectivesBool=True, forceObjectives=False, objectiveFType=None
                                        , optSteps=optStepsRe, optReset=False, pathIsoDoseXML=pathIsoDoseXML, debug=False)
                    print (' - [main] finalPlanStatus: ', finalPlanStatus)
                    if finalPlanStatus:
                        robustEvaluation(planNameFinalAuto, pathRobustTemplate, pathPatient, force=True, verbose=True)
                
                # For evaluation
                if 0:
                    # robustEvaluation(planNameFinal, force=True)
                    # robustEvaluation(planNameFinalAuto, force=True)
                    robustEvaluation(planNameOG, pathRobustTemplate, pathPatient, force=False, verbose=True)
                
                # For NTCP evaluation
                if 0:
                    plansForNTCP = [planNameOG, planNameFinal, planNameFinalAuto]
                    getNTCPVals(patientID, plansForNTCP, pathPatient)

        except:
            traceback.print_exc()
            if DEBUG_PDB: pdb.set_trace()

        print (f' \n\n ===================== end for {pathPatient} [{planNameOG}] (in {round(time.time() - tPatient, 2)} s) ===================== \n\n')

    except:
        traceback.print_exc()
        if DEBUG_PDB: pdb.set_trace()

########################################################
#                       PARAMS(S)                      #
########################################################

if __name__ == "__main__":

    ###################################################################################
    # Step 0 - Init project paths
    ###################################################################################
    DIR_THIS = Path(__file__).parent.absolute() # P:\RayStationScripts
    DIR_DATA = Path(DIR_THIS).parent.absolute().joinpath('RayStationData')
    
    ###################################################################################
    # Step 1 - Logging details
    ###################################################################################
    DIR_LOGS        = Path(DIR_THIS).joinpath('_logs', 'logsProton', 'run4')
    loggerTimestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    pathLogFile     = Path(DIR_LOGS).joinpath("proton-log_{}.txt".format(loggerTimestamp))
    Path(pathLogFile).parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=str(pathLogFile), level=logging.DEBUG, filemode='a', format='%(asctime)s[%(levelname)s] %(funcName)s: %(message)s',datefmt='%d/%m/%Y %I:%M:%S %p')

    ###################################################################################
    # Step 2 - Get common files
    ###################################################################################
    pathKNOProtonsObjectivesClassSolution = Path(DIR_DATA).joinpath('assets', 'objective-template-proton-kno.xml')
    pathRobustEvalTemplate = Path(DIR_DATA).joinpath('assets', 'eval-template-proton-robust.json')
    pathIsoDoseXML = Path(DIR_DATA).joinpath('LUMC-Dose', '_tmp', 'isodose.xml'); # pathIsoDoseXML = None

    ###################################################################################
    # Step 3 -  Specific patient paths (of data extracted from RTPACS)
    ###################################################################################
    # For single patients
    if 0:
        pathPatient                      = Path(DIR_DATA).joinpath('')
        pathKNOProtonsbjectivesClinical = Path(pathPatient).joinpath('.xml')
        
        keynameCancerType = '1P PVPROTONEN'
        optStepsRe        = 4

        contourType = config.KEYNAME_CONTOUR_EVAL # config.{KEYNAME_CONTOUR_EVAL, KEYNAME_CONTOUR_AUTO, KEYNAME_CONTOUR_CLINICAL, KEYNAME_CONTOUR_ALL, KEYNAME_CONTOUR_DEBUG}
    
    # Step 2 - Main
    params = {

        # Patient Case related
        config.KEYNAME_PATH_PATIENT           : pathPatient
        
        # Paths
        , config.KEYNAME_PATH_CLASSSOL        : pathKNOProtonsObjectivesClassSolution
        , config.KEYNAME_PATH_OBJECTIVES      : pathKNOProtonsbjectivesClinical
        , config.KEYNAME_PATH_ROBUST_TEMPLATE : pathRobustEvalTemplate
        , config.KEYNAME_PATH_ISODOSEXML      : pathIsoDoseXML
        
        # Plan parameters
        , config.KEYNAME_CANCER_TYPE            : keynameCancerType
        , config.KEYNAME_OPT_STEPS_RE           : optStepsRe
        , config.KEYNAME_CONTOUR_TYPE           : contourType
    }

    print ('\n -------------------------- [main] params: ')
    print (params)
    main(params)
    
    
# To run code in RS console (and print to console)
"""
import sys
import runpy
from pathlib import Path

sys.path.append(str(Path('H:\\').joinpath('RayStationScripts')))
pathFile = str(Path('H:\\').joinpath('RayStationScripts', 'hnDoseProtons-v1.py'))
_ = runpy.run_path(pathFile, run_name='__main__')

"""