// JavaScript source code
function calculate(exportData) {
    document.body.style.cursor = 'wait';
    var panels = document.getElementsByTagName("*");
    for (i = 0; i < panels.length; i++) {
        panels[i].style.cursor = 'wait';
    }
    setTimeout(function () {
        console.log("Callback Funktion wird aufgerufen");
        calculateNew(exportData);
    }, 10);
}

function calculateNew(exportData) {
    // variable definition

    // general parameters
    var numberSubPopulations = 3;               // number of sub-populations, here three (low risk, high-risk compliant, high-risk non-compliant)
    var zeroPoint;                              // number of days prior to simulation start for which values need to be calculated
    var repeat;
    var capacityNeverExceeded;
    var lastRound;

    // set variations of model parameters
    var betaMod = [];                                       // values for beta used for model uncertainty
    var gammaMod = [];                                      // values for gamma used for model uncertainty

    // variables for each population
    var nHerdImmunity = [];                     // number of immune people in a population necessary for herd immunity
    var herdImmunityReached = [];               // flag indicating whether in a simulation run herd immunity is reached for the population considered
    var herdImmunityAlwaysReached = [];         // flag indicating whether in all simulation runs herd immunity is reached for the population considered
    var iHerdImmunity = [];                     // latest time at which herd immunity is reached in all simulation runs
    var startLevelReached = [];                 // flag indicating whether in a simulation run herd immunity is reached for the population considered and the number of infectious people is below the starting level
    var startLevelAlwaysReached = [];           // flag indicating whether in all simulation runs herd immunity is reached for the population considered and the number of infectious people is below the starting level
    var iStartLevelReached = [];                // latest time at which herd immunity is reached and the number of infectious people is below the starting level in all simulation runs

    // variables for each population and subpopulation
    var subpopulation = [];                     // size of subpopulations of a population (low risk, high-risk compliant, high-risk non-compliant)
    var infRate = [];                           // infection rate of a subpopulations of a population during steady state at the beginning

    // create two dimensional arrays
    for (i = 0; i < numberPopulations; i++) {
        subpopulation[i] = [];
        infRate[i] = [];
    }


    // time dependent variables for all considered sub-populations together
    var totalNewlyNeedingIntensiveCare = [];    // the total number of people newly needing intensive care
    var totalMinNeedingIntensiveCare = [];      // the minimum number of people in the given population in need of intensive care; minimum is meant with respect to the uncertainty of the model
    var totalNeedingIntensiveCare = [];         // the number of people in the given population in need of intensive care
    var totalMaxNeedingIntensiveCare = [];      // the maximum number of people in the given population in need of intensive care; maximum is meant with respect to the uncertainty of the model
    var totalMinHospitalized = [];              // the minimum number of people in the given population that obtain intensive care; minimum is meant with respect to the uncertainty of the model
    var totalHospitalized = [];                 // the number of people in the given population that obtain intensive care
    var totalMaxHospitalized = [];              // the maximum number of people in the given population that obtain intensive care; maximum is meant with respect to the uncertainty of the model
    var hospitalCapacity = [];                  // capacity of intensive care units in the hospitals for the given population

    var totalSusceptible = [];                  // total (i.e. including all sub-populations) number of susceptible individuals in the given population
    var totalInfectious = [];                   // total (i.e. including all sub-populations) number of infectious individuals in the given population
    var totalRemoved = [];                      // total (i.e. including all sub-populations) number of removed (i.e. recovered or deceased) individuals in the given population
    var totalDeltaValue = [];                   // total (i.e. including all sub-populations) number of newly infected individuals in the given population

    var totalMinDeceased = [];                  // the minimum number of deceased people in the given population; minimum is meant with respect to the uncertainty of the model
    var totalDeceased = [];                     // the number of deceased people in the given population
    var totalMaxDeceased = [];                  // the maximum number of deceased people in the given population; maximum is meant with respect to the uncertainty of the model

    var actualContactReduction = [];            // percentage by which the number of contacts is reduced in the given population
    var averageContactReduction = [];           // 14-days average of the percentage by which the number of contacts is reduced in the given population
    var desiredContactReduction = [];           // percentage by which the number of contacts should be reduced according to the policy in place in the given population

    var freeCapacity = [];                      // number of free intensive care units
    var fractionIntensiveCare = [];             // fraction of individuals needing intensive care for which capacity is available

    // create two dimensional arrays
    for (i = 0; i < numberPopulations; i++) {
        totalNewlyNeedingIntensiveCare[i] = [];
        totalMinNeedingIntensiveCare[i] = [];
        totalNeedingIntensiveCare[i] = [];
        totalMaxNeedingIntensiveCare[i] = [];
        totalMinHospitalized[i] = [];
        totalHospitalized[i] = [];
        totalMaxHospitalized[i] = [];
        hospitalCapacity[i] = [];

        totalSusceptible[i] = [];
        totalInfectious[i] = [];
        totalRemoved[i] = [];
        totalDeltaValue[i] = [];

        totalMinDeceased[i] = [];
        totalDeceased[i] = [];
        totalMaxDeceased[i] = [];

        desiredContactReduction[i] = [];
        averageContactReduction[i] = [];
        actualContactReduction[i] = [];

        freeCapacity[i] = [];
        fractionIntensiveCare[i] = [];
    }

    // time dependent variables for subpopulations
    // variables follow the logic: population, subpopulation, time
    var deltaValue = [];                                    // number of newly infected people per population and subpopulation at a given time
    var susceptibles = [];                                  // number of people who have not had the decease yet per population and subpopulation at a given time
    var infectious = [];                                    // number of people who are infected and infectious per population and subpopulation at a given time
    var removed = [];                                       // number of people who are immune or deceased per population and subpopulation at a given time
    var newlyNeedingIntensiveCare = [];                     // number of people who newly need intensive care
    var needingIntensiveCare = [];                          // number of people who need intensive care per population and subpopulation at a given time
    var newlyHospitalized = [];                             // number of people newly admitted to the hospital for intensive care
    var hospitalized = [];                                  // number of people who are in the hospital in intensive care per population and subpopulation at a given time
    var deceased = [];                                      // number of people who have died per population and subpopulation at a given time

    // create two dimensional arrays
    for (i = 0; i < numberPopulations; i++) {
        deltaValue[i] = [];
        susceptibles[i] = [];
        infectious[i] = [];
        removed[i] = [];
        newlyNeedingIntensiveCare[i] = [];
        needingIntensiveCare[i] = [];
        hospitalized[i] = [];
        newlyHospitalized[i] = [];
        deceased[i] = [];

        // create three dimensional arrays
        for (j = 0; j < numberSubPopulations; j++) {
            deltaValue[i][j] = [];
            susceptibles[i][j] = [];
            infectious[i][j] = [];
            removed[i][j] = [];
            newlyNeedingIntensiveCare[i][j] = [];
            needingIntensiveCare[i][j] = [];
            hospitalized[i][j] = [];
            newlyHospitalized[i][j] = [];
            deceased[i][j] = [];
        }
    }

    // variables depending on interaction between populations and subpopulations
    var contactMatrix = [];                                 // matrix defining how the contacts are distributed between the populations and subpopulations
    for (i = 0; i < numberPopulations; i++) {
        contactMatrix[i] = [];
        for (j = 0; j < numberSubPopulations; j++) {
            contactMatrix[i][j] = [];
            for (l = 0; l < numberPopulations; l++) {
                contactMatrix[i][j][l] = [];
            }
        }
    }




    // the time needs to look back at least 14 days (for 14-days average calculation) or the sum of registration delay and duration of hospitalization
    // in order to assure that the start of the simulation is correct
    zeroPoint = Math.max(registrationDelay[0] + durationHospitalization[0], registrationDelay[1] + durationHospitalization[1], 14);

    // initialize values
    // initialize all values to 0 exect the "min" values which are set to the maximum possible value
    for (j = 0; j < numberPopulations; j++) {
        for (i = 0; i < durationSimulation[j] + 2 * zeroPoint; i++) {

            for (k = 0; k < numberSubPopulations; k++) {
                deltaValue[j][k][i] = 0;
                susceptibles[j][k][i] = 0;
                infectious[j][k][i] = 0;
                removed[j][k][i] = 0;
                newlyNeedingIntensiveCare[j][k][i] = 0;
                needingIntensiveCare[j][k][i] = 0;
                hospitalized[j][k][i] = 0;
                newlyHospitalized[j][k][i] = 0;
                deceased[j][k][i] = 0;
            }

            totalMinHospitalized[j][i] = population[j];
            totalHospitalized[j][i] = 0;
            totalMaxHospitalized[j][i] = 0;

            totalNewlyNeedingIntensiveCare[j][i] = 0;
            totalMinNeedingIntensiveCare[j][i] = population[j];
            totalNeedingIntensiveCare[j][i] = 0;
            totalMaxNeedingIntensiveCare[j][i] = 0;

            totalSusceptible[j][i] = 0;
            totalInfectious[j][i] = 0;
            totalRemoved[j][i] = 0;

            desiredContactReduction[j][i] = 0;
            averageContactReduction[j][i] = 0;
            actualContactReduction[j][i] = 0;

            totalMinDeceased[j][i] = population[j];
            totalDeceased[j][i] = 0;
            totalMaxDeceased[j][i] = 0;
            totalDeltaValue[j][i] = 0;
        }
    }

    // determine initial values, i.e. data at the zero point: translate available data into initial simulation variables
    for (i = 0; i < numberPopulations; i++) {

        for (j = 0; j < numberSubPopulations; j++) {
            // low risk
            if (j == 0) {
                deceased[i][j][zeroPoint] = initiallyRegisteredDeceased[i] * (100 - percentageHighRisk[i]) * riskIntensiveCareLowRisk[i] / (riskIntensiveCareLowRisk[i] * (100 - percentageHighRisk[i]) + riskIntensiveCareHighRisk[i] * percentageHighRisk[i]);
                removed[i][j][zeroPoint] = initiallyRegisteredRecovered[i] * factorUndercoverage[i] * (100 - percentageHighRisk[i]) / 100 + deceased[i][j][zeroPoint];
                infectious[i][j][zeroPoint] = ((initiallyRegisteredCases[i] - initiallyRegisteredRecovered[i]) * factorUndercoverage[i] - initiallyRegisteredDeceased[i]) * (100 - percentageHighRisk[i]) / 100;
                subpopulation[i][j] = population[i] * (100 - percentageHighRisk[i]) / 100;
            }
            // high risk compliant
            if (j == 1) {
                deceased[i][j][zeroPoint] = initiallyRegisteredDeceased[i] * (percentageHighRisk[i] * complianceHighRiskGroup[i]) / 100 * riskIntensiveCareHighRisk[i] / (riskIntensiveCareLowRisk[i] * (100 - percentageHighRisk[i]) + riskIntensiveCareHighRisk[i] * percentageHighRisk[i]);
                removed[i][j][zeroPoint] = initiallyRegisteredRecovered[i] * factorUndercoverage[i] * percentageHighRisk[i] * complianceHighRiskGroup[i] / (100 * 100) + deceased[i][j][zeroPoint];
                infectious[i][j][zeroPoint] = ((initiallyRegisteredCases[i] - initiallyRegisteredRecovered[i]) * factorUndercoverage[i] - initiallyRegisteredDeceased[i]) * percentageHighRisk[i] * complianceHighRiskGroup[i] / (100 * 100);
                subpopulation[i][j] = population[i] * percentageHighRisk[i] / 100 * complianceHighRiskGroup[i] / 100;
            }
            // high risk non-compliant
            if (j == 2) {
                deceased[i][j][zeroPoint] = initiallyRegisteredDeceased[i] * (percentageHighRisk[i] * (100 - complianceHighRiskGroup[i])) / 100 * riskIntensiveCareHighRisk[i] / (riskIntensiveCareLowRisk[i] * (100 - percentageHighRisk[i]) + riskIntensiveCareHighRisk[i] * percentageHighRisk[i]);
                removed[i][j][zeroPoint] = initiallyRegisteredRecovered[i] * factorUndercoverage[i] * (percentageHighRisk[i] * (100 - complianceHighRiskGroup[i])) / (100 * 100) + deceased[i][j][zeroPoint];
                infectious[i][j][zeroPoint] = ((initiallyRegisteredCases[i] - initiallyRegisteredRecovered[i]) * factorUndercoverage[i] - initiallyRegisteredDeceased[i]) * percentageHighRisk[i] * (100 - complianceHighRiskGroup[i]) / (100 * 100);
                subpopulation[i][j] = population[i] * percentageHighRisk[i] / 100 * (100 - complianceHighRiskGroup[i]) / 100;
            }

            susceptibles[i][j][zeroPoint] = subpopulation[i][j] - infectious[i][j][zeroPoint] - removed[i][j][zeroPoint];

            // check for values below 0
            if (susceptibles[i][j][zeroPoint] < 0) susceptibles[i][j][zeroPoint] = 0;

        }
        // we start the simulation out of a lockdown situation, i.e. we set all contact levels to low
        actualContactReduction[i][zeroPoint] = (100.0 - lowContactLevel[i]);
        desiredContactReduction[i][zeroPoint] = (100.0 - lowContactLevel[i]);
        averageContactReduction[i][zeroPoint] = (100.0 - lowContactLevel[i]);
    }

    // determine capacity of health care system for all populations
    // capacity is constant up to zero point, then increases as specified
    for (j = 0; j < numberPopulations; j++) {
        for (i = 0; i < durationSimulation[j] + 2 * zeroPoint; i++) {
            if (i > zeroPoint)
                hospitalCapacity[j][i] = intensiveCareUnits[j] + dailyIncreaseIntensiveCareUnits[j] * (i - zeroPoint);
            else
                hospitalCapacity[j][i] = intensiveCareUnits[j]
        }
    }

    // initialize values for herd immunity and infection risk
    for (j = 0; j < numberPopulations; j++) {
        herdImmunityAlwaysReached[j] = true;
        startLevelAlwaysReached[j] = true;
        iHerdImmunity[j] = 0;
        iStartLevelReached[j] = 0;
    }


    // run simulation

    // calculate backwards for registrationsDelay timesteps
    // assume constant infection rate

    // determine infection rate
    for (i = 0; i < numberPopulations; i++) {
        totalNeedingIntensiveCare[i][zeroPoint] = 0;
        for (j = 0; j < numberSubPopulations; j++) {
            // during steady state, the infection rate is equal to the recovery rate, which is easy to calculate
            infRate[i][j] = gamma[i] / 100.0 * infectious[i][j][zeroPoint];
            // at steady state, the number of requires intensive care units equals the infection rate times the risk for needing intensive care time number of days of a hospital stay. The risk needs to be adjusted for the factor for undercoverage
            // we start by assuming that we have sufficient capacity in the hospitals. If not, we will correct
            if (j == 0) {
                newlyNeedingIntensiveCare[i][j][zeroPoint] = infRate[i][j] * riskIntensiveCareLowRisk[i] / 100.0 / factorUndercoverage[i];
                needingIntensiveCare[i][j][zeroPoint] = durationHospitalization[i] * infRate[i][j] * riskIntensiveCareLowRisk[i] / 100.0 / factorUndercoverage[i];
                hospitalized[i][j][zeroPoint] = needingIntensiveCare[i][j][zeroPoint];
                newlyHospitalized[i][j][zeroPoint] = newlyNeedingIntensiveCare[i][j][zeroPoint];
                totalHospitalized[i][zeroPoint] += hospitalized[i][j][zeroPoint];
            }
            else {
                newlyNeedingIntensiveCare[i][j][zeroPoint] = infRate[i][j] * riskIntensiveCareHighRisk[i] / 100.0 / factorUndercoverage[i];
                needingIntensiveCare[i][j][zeroPoint] = durationHospitalization[i] * infRate[i][j] * riskIntensiveCareHighRisk[i] / 100.0 / factorUndercoverage[i];
                hospitalized[i][j][zeroPoint] = needingIntensiveCare[i][j][zeroPoint];
                newlyHospitalized[i][j][zeroPoint] = newlyNeedingIntensiveCare[i][j][zeroPoint];
                totalHospitalized[i][zeroPoint] += hospitalized[i][j][zeroPoint];
            }
            totalNewlyNeedingIntensiveCare[i][zeroPoint] += newlyNeedingIntensiveCare[i][j][zeroPoint];
            totalNeedingIntensiveCare[i][zeroPoint] += needingIntensiveCare[i][j][zeroPoint];
        }
    }

    // if the total number of people needing intensive care exceeds the capacity, we assume that the hospitaized patients are representative of the population, i.e. no preference for any group when assigning hospital beds
    for (i = 0; i < numberPopulations; i++) {
        if (totalNeedingIntensiveCare[i][zeroPoint] > hospitalCapacity[i][zeroPoint]) {
            totalHospitalized[i][zeroPoint] = 0;
            for (j = 0; j < numberSubPopulations; j++) {
                newlyHospitalized[i][j][zeroPoint] = hospitalCapacity[i][zeroPoint] * subpopulation[i][j] / population[i] / durationHospitalization[i];
                hospitalized[i][j][zeroPoint] = hospitalCapacity[i][zeroPoint] * subpopulation[i][j] / population[i];
                totalHospitalized[i][zeroPoint] += hospitalized[i][j][zeroPoint];
            }
        }
    }

    // determine free capacity and fractionof free capacity
    for (j = 0; j < numberPopulations; j++) {
        freeCapacity[j][zeroPoint] = Math.max(hospitalCapacity[j] - totalHospitalized[j][zeroPoint], 0);
        // if nobody needs intensive care, 1 = 100% of the capacity is free, i.e. 100% of the patients needing intesive care can be given intensive care
        if (totalNeedingIntensiveCare[j][zeroPoint] == 0) fractionIntensiveCare[j][zeroPoint] = 1;
        // if more people need intensive care than capacity is free, only a fraction can be given intensive care
        else fractionIntensiveCare[j][zeroPoint] = Math.min(freeCapacity[j][zeroPoint] / totalNeedingIntensiveCare[j][zeroPoint], 1);
    }

    // calculate backward using steady state assumption
    for (i = zeroPoint - 1; i >= 0; i--) {
        for (j = 0; j < numberPopulations; j++) {

            // we start out of lockdown
            actualContactReduction[j][i] = (100.0 - lowContactLevel[j]);
            desiredContactReduction[j][i] = (100.0 - lowContactLevel[j]);
            averageContactReduction[j][i] = (100.0 - lowContactLevel[j]);

            // calculate values
            for (k = 0; k < numberSubPopulations; k++) {
                deltaValue[j][k][i] = infRate[j][k];
                susceptibles[j][k][i] = susceptibles[j][k][i + 1] + infRate[j][k];
                infectious[j][k][i] = infectious[j][k][i + 1];
                removed[j][k][i] = removed[j][k][i + 1] - infRate[j][k];
                newlyNeedingIntensiveCare[j][k][i] = newlyNeedingIntensiveCare[j][k][i + 1];
                needingIntensiveCare[j][k][i] = needingIntensiveCare[j][k][i + 1];
                hospitalized[j][k][i] = hospitalized[j][k][i + 1];
                newlyHospitalized[j][k][i] = newlyHospitalized[j][k][i + 1];
            }

            // calculate total number of people needing intensive care in the given population at the given time
            for (k = 0; k < numberSubPopulations; k++) {
                totalNeedingIntensiveCare[j][i] += needingIntensiveCare[j][k][i];
            }

            // determine free capacity and fraction of free capacity
            freeCapacity[j][i] = Math.max(hospitalCapacity[j][i] - totalNeedingIntensiveCare[j][i], 0);
            // if nobody needs intensive care, 1 = 100% of the capacity is free, i.e. 100% of the patients needing intesive care can be given intensive care
            if (totalNeedingIntensiveCare[j][i] == 0) fractionIntensiveCare[j][i] = 1;
            // if more people need intensive care than capacity is free, only a fraction can be given intensive care
            else fractionIntensiveCare[j][i] = Math.min(hospitalCapacity[j][i] / totalNeedingIntensiveCare[j][i], 1);


            for (k = 0; k < numberSubPopulations; k++) {
                deceased[j][k][i] = deceased[j][k][i + 1];
                if (k == 0) {
                    // the ones who need intensive care for which intensive care is available die with the probability defined, the others die for sure
                    // the risk needs to be adjusted for the factor for undercoverage
                    deceased[j][k][i] -= infRate[j][k] * riskIntensiveCareLowRisk[j] / 100.0 * caseMortality[j] / 100.0 * fractionIntensiveCare[j][i] / factorUndercoverage[j];
                    deceased[j][k][i] -= infRate[j][k] * riskIntensiveCareLowRisk[j] / 100.0 * (1 - fractionIntensiveCare[j][i]) / factorUndercoverage[j];
                }
                else {
                    deceased[j][k][i] -= infRate[j][k] * riskIntensiveCareHighRisk[j] / 100.0 * caseMortality[j] / 100.0 * fractionIntensiveCare[j][i] / factorUndercoverage[j];
                    deceased[j][k][i] -= infRate[j][k] * riskIntensiveCareHighRisk[j] / 100.0 * (1 - fractionIntensiveCare[j][i]) / factorUndercoverage[j];
                }
            }
        }
    }

    // calculate aggregate values for display
    for (i = zeroPoint; i >= 0; i--) {
        for (j = 0; j < numberPopulations; j++) {
            totalSusceptible[j][i] = 0;
            totalInfectious[j][i] = 0;
            totalRemoved[j][i] = 0;
            totalHospitalized[j][i] = 0;
            totalDeceased[j][i] = 0;

            for (k = 0; k < numberSubPopulations; k++) {
                deltaValue[j][k][i] = infRate[j][k];
                totalSusceptible[j][i] += susceptibles[j][k][i];
                totalInfectious[j][i] += infectious[j][k][i];
                totalRemoved[j][i] += removed[j][k][i];
                totalDeceased[j][i] += deceased[j][k][i];
                totalHospitalized[j][i] += hospitalized[j][k][i];

            }
            // determine the min and max values for display
            if (totalDeceased[j][i] < totalMinDeceased[j][i]) totalMinDeceased[j][i] = totalDeceased[j][i];
            if (totalDeceased[j][i] > totalMaxDeceased[j][i]) totalMaxDeceased[j][i] = totalDeceased[j][i];
            if (totalHospitalized[j][i] < totalMinHospitalized[j][i]) totalMinHospitalized[j][i] = totalHospitalized[j][i];
            if (totalHospitalized[j][i] > totalMaxHospitalized[j][i]) totalMaxHospitalized[j][i] = totalHospitalized[j][i];
            if (totalNeedingIntensiveCare[j][i] < totalMinNeedingIntensiveCare[j][i]) totalMinNeedingIntensiveCare[j][i] = totalNeedingIntensiveCare[j][i];
            if (totalNeedingIntensiveCare[j][i] > totalMaxNeedingIntensiveCare[j][i]) totalMaxNeedingIntensiveCare[j][i] = totalNeedingIntensiveCare[j][i];
        }
    }

    repeat = true;
    capacityNeverExceeded = true;
    lastRound = true;

    if (strategy[populationIndex] == 2) {
        highContactLevel[populationIndex] = 0;
        lowContactLevel[populationIndex] = 0;
        lastRound = false;
    }

    console.log(numberSimulations[populationIndex]);

    // forward calculation
    while (repeat) {
        if (lastRound) repeat = false;
        for (j = 0; j < numberPopulations; j++) {
            for (i = zeroPoint + 1; i < durationSimulation[0] + 2 * zeroPoint; i++) {

                totalMinDeceased[j][i] = population[j];
                totalDeceased[j][i] = 0;
                totalMaxDeceased[j][i] = 0;
                totalMinHospitalized[j][i] = population[j];
                totalHospitalized[j][i]
                totalMaxHospitalized[j][i] = 0;
                totalMinNeedingIntensiveCare[j][i] = population[j];
                totalNeedingIntensiveCare[j][i] = 0;
                totalMaxNeedingIntensiveCare[j][i] = 0;
            }
        }
        for (round = 0; round < numberSimulations[populationIndex]; round++) {


            // calculate for all considered model variants
            for (modelIndex = 0; modelIndex < 5; modelIndex++) {

                // initialize model specific variables for herd immunity and infection level
                for (j = 0; j < numberPopulations; j++) {
                    herdImmunityReached[j] = false;
                    startLevelReached[j] = false;
                }

                // loop over all future days
                for (i = zeroPoint + 1; i < durationSimulation[0] + zeroPoint + 1; i++) {

                    // loop over all populations
                    for (j = 0; j < numberPopulations; j++) {


                        // at the beginning, calculate the contact matrix based on parameters and options
                        if (i == zeroPoint + 1) calculateContactMatrix(j, option[j]);



                        // if interventions are over, calculate contact matrix based on parameters without options for special care for the high risk group
                        if (i == durationIntervention[j] + zeroPoint + 1) {
                            calculateContactMatrix(j, 0);
                        }

                        // choose model parameters to represent uncertainty of model
                        betaMod[0] = parseFloat(beta[j] * (1 + uncertaintyModel[j] / 100.0));
                        betaMod[1] = parseFloat(beta[j] * (1 - uncertaintyModel[j] / 100.0));
                        betaMod[2] = parseFloat(beta[j]);
                        betaMod[3] = parseFloat(beta[j]);
                        betaMod[4] = parseFloat(beta[j]);

                        gammaMod[0] = parseFloat(gamma[j]);
                        gammaMod[1] = parseFloat(gamma[j]);
                        gammaMod[2] = parseFloat(gamma[j] * (1 + uncertaintyModel[j] / 100.0));
                        gammaMod[3] = parseFloat(gamma[j] * (1 - uncertaintyModel[j] / 100.0));
                        gammaMod[4] = parseFloat(gamma[j]);

                        // calculate herd immunity number based on model and chosen parameters
                        nHerdImmunity[j] = population[j] - gammaMod[modelIndex] * population[j] / (betaMod[modelIndex] * contactLevelAfterIntervention[j]);

                        // set aggegate variables to 0
                        totalSusceptible[j][i] = 0;
                        totalInfectious[j][i] = 0;
                        totalRemoved[j][i] = 0;
                        totalNeedingIntensiveCare[j][i] = 0;
                        totalNewlyNeedingIntensiveCare[j][i] = 0;
                        totalHospitalized[j][i] = 0;
                        totalDeltaValue[j][i] = 0;

                        // determine contact reduction

                        // strategy 0: no action
                        if (strategy[j] == 0) desiredContactReduction[j][i] = 0;
                        // strategy 1: mitigation
                        if (strategy[j] == 1) desiredContactReduction[j][i] = 100 - lowContactLevel[j];
                        // strategy 2: flatten the curve
                        if (strategy[j] == 2) {
                            desiredContactReduction[j][i] = 100 - lowContactLevel[j];
                        }
                        // strategy 3: adaptive triggering
                        if (strategy[j] == 3) {
                            if (i < 7 + registrationDelay[j]) desiredContactReduction[j][i] = 100 - lowContactLevel[j];
                            else {
                                // count new infections in last 7 days
                                newInfections = 0;
                                for (l = 0; l < 7; l++) newInfections += totalDeltaValue[j][i - l - registrationDelay[j]];
                                if (newInfections / factorUndercoverage[j] < (population[j] / 100000) * thresholdInfections[j]) desiredContactReduction[j][i] = 100 - highContactLevel[j];
                                else desiredContactReduction[j][i] = 100 - lowContactLevel[j];
                            }
                        }
                        // strategy 4: maximum contact while respecting limit of health care system
                        // calculate as many days intot the future as it takes for infected people to present themselves to the hospital
                        // use the past decisions and evaluate both options (high or low contact level) on the given day
                        // do this for all sets of parameters assuming the worst casse contact level
                        if (strategy[j] == 4) {

                            capacityExceeded = false;

                            for (n = 0; n < 5; n++) {
                                for (p = 0; p < Math.min(registrationDelay[j] + durationHospitalization[j], durationSimulation[j] + zeroPoint - i); p++) {
                                    totalNewlyNeedingIntensiveCare[j][i + p] = 0;


                                    // loop over all sub-populations
                                    for (k = 0; k < numberSubPopulations; k++) {

                                        if (p == 0) {
                                            actualContactReduction[j][i + p] = 100 - highContactLevel[j] - uncertaintyContactLevel[j];
                                        }
                                        else {
                                            actualContactReduction[j][i + p] = 100 - lowContactLevel[j] - uncertaintyContactLevel[j];
                                        }

                                        if (actualContactReduction[j][i + p] < 0) actualContactReduction[j][i + p] = 0;
                                        if (actualContactReduction[j][i + p] > 100) actualContactReduction[j][i + p] = 100;

                                        susceptibles[j][k][i + p] = susceptibles[j][k][i - 1 + p];
                                        infectious[j][k][i + p] = infectious[j][k][i - 1 + p];
                                        removed[j][k][i + p] = removed[j][k][i - 1 + p];

                                        // initialize newly infected individuals to 0
                                        deltaValue[j][k][i + p] = 0;

                                        // calculate newly infected individuals through transmission from all populations and subpopulations
                                        for (l = 0; l < numberPopulations; l++) {
                                            for (m = 0; m < numberSubPopulations; m++) {
                                                if (subpopulation[l][m] != 0) {
                                                    deltaValue[j][k][i + p] += (seasonalityMinPercentage[j] + (100 - seasonalityMinPercentage[j]) / 2 + (100 - seasonalityMinPercentage[j]) / 2 * Math.cos((i + p - registrationDelay[j] - zeroPoint - seasonalityPeakDate[j]) / 365 * 2 * Math.PI)) / 100.0
                                                        * betaMod[n] * susceptibles[j][k][i - 1 + p] * contactMatrix[j][k][l][m] * (100 - actualContactReduction[j][i + p]) / 100.0 * infectious[l][m][i - 1 + p] / subpopulation[l][m];
                                                }
                                            }
                                        }

                                        // update SIR model
                                        susceptibles[j][k][i + p] -= deltaValue[j][k][i + p];
                                        infectious[j][k][i + p] += deltaValue[j][k][i + p];
                                        infectious[j][k][i + p] -= gammaMod[n] / 100.0 * infectious[j][k][i - 1 + p];
                                        removed[j][k][i + p] += gammaMod[n] / 100.0 * infectious[j][k][i - 1 + p];

                                        // check for plausibility (number of infectious individuals cannot be lower than 0)
                                        if (infectious[j][k][i + p] < 0) infectious[j][k][i + p] = 0;

                                        // calculate aggregates over sub-populations
                                        totalSusceptible[j][i + p] += susceptibles[j][k][i + p];
                                        totalInfectious[j][i + p] += infectious[j][k][i + p];
                                        totalRemoved[j][i + p] += removed[j][k][i + p];
                                        totalDeltaValue[j][i + p] += deltaValue[j][k][i + p];

                                    }

                                    // initialize values for hospitalized and needing intensive care to previous values
                                    for (k = 0; k < numberSubPopulations; k++) {
                                        hospitalized[j][k][i + p] = hospitalized[j][k][i - 1 + p];
                                        needingIntensiveCare[j][k][i + p] = needingIntensiveCare[j][k][i - 1 + p];
                                    }

                                    // determine number of patients newly needing intensive care and add those to the number of patients needing intensive care
                                    for (k = 0; k < numberSubPopulations; k++) {
                                        if (k == 0) {
                                            newlyNeedingIntensiveCare[j][k][i + p] = deltaValue[j][k][i + p - registrationDelay[j]] * riskIntensiveCareLowRisk[j] / 100.0 / factorUndercoverage[j];
                                        }
                                        if (k > 0) {
                                            newlyNeedingIntensiveCare[j][k][i + p] = deltaValue[j][k][i + p - registrationDelay[j]] * riskIntensiveCareHighRisk[j] / 100.0 / factorUndercoverage[j];
                                        }
                                        totalNewlyNeedingIntensiveCare[j][i + p] += newlyNeedingIntensiveCare[j][k][i + p];
                                        needingIntensiveCare[j][k][i + p] += newlyNeedingIntensiveCare[j][k][i + p];
                                    }

                                    // remove from hospitalized the patients that were admitted at the current time reduced by the defined average length of a hospital stay
                                    for (k = 0; k < numberSubPopulations; k++) {
                                        hospitalized[j][k][i + p] -= newlyHospitalized[j][k][i + p - durationHospitalization[j]];
                                        needingIntensiveCare[j][k][i + p] -= newlyNeedingIntensiveCare[j][k][i + p - durationHospitalization[j]];
                                    }

                                    // determine total number of hospitalized and needing intensive care
                                    totalHospitalized[j][i + p] = 0;
                                    totalNeedingIntensiveCare[j][i + p] = 0;
                                    for (k = 0; k < numberSubPopulations; k++) {
                                        totalHospitalized[j][i + p] += hospitalized[j][k][i + p];
                                        totalNeedingIntensiveCare[j][i + p] += needingIntensiveCare[j][k][i + p];
                                    }

                                    // determine free capacity
                                    freeCapacity[j][i + p] = Math.max(hospitalCapacity[j][i + p - registrationDelay[j]] - totalHospitalized[j][i + p], 0);

                                    // determine fraction of patients that can obtain intensive care
                                    if (totalNeedingIntensiveCare[j][i + p] == 0) fractionIntensiveCare[j][i + p] = 1;
                                    else fractionIntensiveCare[j][i + p] = Math.min(freeCapacity[j][i + p] / totalNewlyNeedingIntensiveCare[j][i + p], 1);



                                    // determine what happens to the newly needing intensive care: admission to intensive care or will they die
                                    for (k = 0; k < numberSubPopulations; k++) {

                                        deceased[j][k][i + p] = deceased[j][k][i - 1 + p];

                                        // use fraction of free capacity to assign mortality rate
                                        if (k == 0) deceased[j][k][i + p] += deltaValue[j][k][i + p - registrationDelay[j]] * riskIntensiveCareLowRisk[j] / 100.0 * caseMortality[j] / 100.0 * fractionIntensiveCare[j][i + p] / factorUndercoverage[j];
                                        if (k == 0) deceased[j][k][i + p] += deltaValue[j][k][i + p - registrationDelay[j]] * riskIntensiveCareLowRisk[j] / 100.0 * (1 - fractionIntensiveCare[j][i + p]) / factorUndercoverage[j];
                                        if (k > 0) deceased[j][k][i + p] += deltaValue[j][k][i + p - registrationDelay[j]] * riskIntensiveCareHighRisk[j] / 100.0 * caseMortality[j] / 100.0 * fractionIntensiveCare[j][i + p] / factorUndercoverage[j];
                                        if (k > 0) deceased[j][k][i + p] += deltaValue[j][k][i + p - registrationDelay[j]] * riskIntensiveCareHighRisk[j] / 100.0 * (1 - fractionIntensiveCare[j][i + p]) / factorUndercoverage[j];

                                        newlyHospitalized[j][k][i + p] = newlyNeedingIntensiveCare[j][k][i + p] * fractionIntensiveCare[j][i + p];
                                        hospitalized[j][k][i + p] += newlyHospitalized[j][k][i + p];
                                    }

                                    totalHospitalized[j][i + p] = 0;
                                    totalDeceased[j][i + p] = 0;
                                    // determine total numer of hospitalized and deceased for display
                                    for (k = 0; k < numberSubPopulations; k++) {
                                        totalHospitalized[j][i + p] += hospitalized[j][k][i + p];
                                        totalDeceased[j][i + p] += deceased[j][k][i + p];
                                    }
                                    // check if there is enough place in intensive care
                                    if (fractionIntensiveCare[j][i + p] < 1) {
                                        capacityExceeded = true;

                                        p = Math.min(registrationDelay[j] + durationHospitalization[j], durationSimulation[j] + zeroPoint - i);
                                    }
                                }
                            }
                            if (capacityExceeded) desiredContactReduction[j][i] = 100 - lowContactLevel[j];
                            else desiredContactReduction[j][i] = 100 - highContactLevel[j];
                        }

                        // set aggegate variables to 0
                        totalSusceptible[j][i] = 0;
                        totalInfectious[j][i] = 0;
                        totalRemoved[j][i] = 0;
                        totalNeedingIntensiveCare[j][i] = 0;
                        totalNewlyNeedingIntensiveCare[j][i] = 0;
                        totalHospitalized[j][i] = 0;
                        totalDeltaValue[j][i] = 0;

                        // if interventation is over, revert to final contact level
                        if (i > durationIntervention[j] + zeroPoint) {
                            desiredContactReduction[j][i] = 100 - contactLevelAfterIntervention[j];
                        }

                        // add noise to desired contact reduction
                        actualContactReduction[j][i] = desiredContactReduction[j][i] + (Math.random() - 0.5) * uncertaintyContactLevel[j];

                        // assure that contact reduction is not < 0 and not > 100
                        if (actualContactReduction[j][i] < 0) actualContactReduction[j][i] = 0;
                        if (actualContactReduction[j][i] > 100) actualContactReduction[j][i] = 100;

                        // calculate average contact reduction
                        averageContactReduction[j][i] = 0;
                        for (l = 0; l < 14; l++)
                            averageContactReduction[j][i] += actualContactReduction[j][i - l] / 14.0;

                        // loop over all sub-populations
                        for (k = 0; k < numberSubPopulations; k++) {

                            susceptibles[j][k][i] = susceptibles[j][k][i - 1];
                            infectious[j][k][i] = infectious[j][k][i - 1];
                            removed[j][k][i] = removed[j][k][i - 1];

                            // initialize newly infected individuals to 0
                            deltaValue[j][k][i] = 0;

                            // calculate newly infected individuals through transmission from all populations and subpopulations
                            for (l = 0; l < numberPopulations; l++) {
                                for (m = 0; m < numberSubPopulations; m++) {
                                    if (subpopulation[l][m] != 0) {
                                        deltaValue[j][k][i] += (seasonalityMinPercentage[j] + (100 - seasonalityMinPercentage[j]) / 2 + (100 - seasonalityMinPercentage[j]) / 2 * Math.cos((i - zeroPoint - seasonalityPeakDate[j]) / 365 * 2 * Math.PI)) / 100.0
                                            * betaMod[modelIndex] * susceptibles[j][k][i - 1] * contactMatrix[j][k][l][m] * (100 - actualContactReduction[j][i]) / 100.0 * infectious[l][m][i - 1] / subpopulation[l][m];
                                    }
                                }
                            }

                            // update SIR model
                            susceptibles[j][k][i] -= deltaValue[j][k][i];
                            infectious[j][k][i] += deltaValue[j][k][i];
                            infectious[j][k][i] -= gammaMod[modelIndex] / 100.0 * infectious[j][k][i - 1];
                            removed[j][k][i] += gammaMod[modelIndex] / 100.0 * infectious[j][k][i - 1];

                            // check for plausibility (number of infectious individuals cannot be lower than 0)
                            if (infectious[j][k][i] < 0) infectious[j][k][i] = 0;

                            // calculate aggregates over sub-populations
                            totalSusceptible[j][i] += susceptibles[j][k][i];
                            totalInfectious[j][i] += infectious[j][k][i];
                            totalRemoved[j][i] += removed[j][k][i];
                            totalDeltaValue[j][i] += deltaValue[j][k][i];

                        }

                        // initialize values for hospitalized and needing intensive care to previous values
                        for (k = 0; k < numberSubPopulations; k++) {
                            hospitalized[j][k][i] = hospitalized[j][k][i - 1];
                            needingIntensiveCare[j][k][i] = needingIntensiveCare[j][k][i - 1];
                        }

                        // determine number of patients newly needing intensive care and add those to the number of patients needing intensive care
                        for (k = 0; k < numberSubPopulations; k++) {
                            if (k == 0) {
                                newlyNeedingIntensiveCare[j][k][i] = deltaValue[j][k][i - registrationDelay[j]] * riskIntensiveCareLowRisk[j] / 100.0 / factorUndercoverage[j];
                            }
                            if (k > 0) {
                                newlyNeedingIntensiveCare[j][k][i] = deltaValue[j][k][i - registrationDelay[j]] * riskIntensiveCareHighRisk[j] / 100.0 / factorUndercoverage[j];
                            }
                            totalNewlyNeedingIntensiveCare[j][i] += newlyNeedingIntensiveCare[j][k][i];
                            needingIntensiveCare[j][k][i] += newlyNeedingIntensiveCare[j][k][i];
                        }

                        // remove from hospitalized the patients that were admitted at the current time reduced by the defined average length of a hospital stay
                        for (k = 0; k < numberSubPopulations; k++) {
                            hospitalized[j][k][i] -= newlyHospitalized[j][k][i - durationHospitalization[j]];
                            needingIntensiveCare[j][k][i] -= newlyNeedingIntensiveCare[j][k][i - durationHospitalization[j]];
                        }

                        // determine total number of hospitalized and needing intensive care
                        totalHospitalized[j][i] = 0;
                        totalNeedingIntensiveCare[j][i] = 0;
                        for (k = 0; k < numberSubPopulations; k++) {
                            totalHospitalized[j][i] += hospitalized[j][k][i];
                            totalNeedingIntensiveCare[j][i] += needingIntensiveCare[j][k][i];
                        }

                        // determine free capacity
                        freeCapacity[j][i] = Math.max(hospitalCapacity[j][i] - totalHospitalized[j][i], 0);

                        // determine fraction of patients that can obtain intensive care
                        if (totalNeedingIntensiveCare[j][i] == 0) fractionIntensiveCare[j][i] = 1;
                        else fractionIntensiveCare[j][i] = Math.min(freeCapacity[j][i] / totalNewlyNeedingIntensiveCare[j][i], 1);

                        if (i < durationIntervention[populationIndex] + zeroPoint && fractionIntensiveCare[populationIndex][i] < 1) capacityNeverExceeded = false;

                        // determine what happens to the newly needing intensive care: admission to intensive care or will they die
                        for (k = 0; k < numberSubPopulations; k++) {
                            deceased[j][k][i] = deceased[j][k][i - 1];
                            // use fraction of free capacity to assign mortality rate
                            if (k == 0) deceased[j][k][i] += deltaValue[j][k][i - registrationDelay[j]] * riskIntensiveCareLowRisk[j] / 100.0 * caseMortality[j] / 100.0 * fractionIntensiveCare[j][i] / factorUndercoverage[j];
                            if (k == 0) deceased[j][k][i] += deltaValue[j][k][i - registrationDelay[j]] * riskIntensiveCareLowRisk[j] / 100.0 * (1 - fractionIntensiveCare[j][i]) / factorUndercoverage[j];
                            if (k > 0) deceased[j][k][i] += deltaValue[j][k][i - registrationDelay[j]] * riskIntensiveCareHighRisk[j] / 100.0 * caseMortality[j] / 100.0 * fractionIntensiveCare[j][i] / factorUndercoverage[j];
                            if (k > 0) deceased[j][k][i] += deltaValue[j][k][i - registrationDelay[j]] * riskIntensiveCareHighRisk[j] / 100.0 * (1 - fractionIntensiveCare[j][i]) / factorUndercoverage[j];

                            newlyHospitalized[j][k][i] = newlyNeedingIntensiveCare[j][k][i] * fractionIntensiveCare[j][i];
                            hospitalized[j][k][i] += newlyHospitalized[j][k][i];
                        }

                        totalHospitalized[j][i] = 0;
                        totalDeceased[j][i] = 0;

                        // determine total numer of hospitalized and deceased for display
                        for (k = 0; k < numberSubPopulations; k++) {
                            totalHospitalized[j][i] += hospitalized[j][k][i];
                            totalDeceased[j][i] += deceased[j][k][i];
                        }

                        // determine min and max of totals for display
                        if (totalDeceased[j][i] < totalMinDeceased[j][i]) totalMinDeceased[j][i] = totalDeceased[j][i];
                        if (totalDeceased[j][i] > totalMaxDeceased[j][i]) totalMaxDeceased[j][i] = totalDeceased[j][i];
                        if (totalHospitalized[j][i] < totalMinHospitalized[j][i]) totalMinHospitalized[j][i] = totalHospitalized[j][i];
                        if (totalHospitalized[j][i] > totalMaxHospitalized[j][i]) totalMaxHospitalized[j][i] = totalHospitalized[j][i];
                        if (totalNeedingIntensiveCare[j][i] < totalMinNeedingIntensiveCare[j][i]) totalMinNeedingIntensiveCare[j][i] = totalNeedingIntensiveCare[j][i];
                        if (totalNeedingIntensiveCare[j][i] > totalMaxNeedingIntensiveCare[j][i]) totalMaxNeedingIntensiveCare[j][i] = totalNeedingIntensiveCare[j][i];

                        // determine events for display

                        // determine whether herd immunity is reached
                        if (population[j] - totalSusceptible[j][i] > nHerdImmunity[j] && !herdImmunityReached[j]) {
                            herdImmunityReached[j] = true;
                            if (i > iHerdImmunity[j]) iHerdImmunity[j] = i;
                        }

                        // determine whether initial risk of infection is reached
                        if (totalInfectious[j][i] < totalInfectious[j][zeroPoint] && herdImmunityReached[j] && !startLevelReached[j]) {
                            startLevelReached[j] = true;
                            if (i > iStartLevelReached[j]) iStartLevelReached[j] = i;
                        }
                    }
                }
            }
            for (j = 0; j < numberPopulations; j++) {
                if (!herdImmunityReached[j]) herdImmunityAlwaysReached[j] = false;
                if (!startLevelReached[j]) startLevelAlwaysReached[j] = false;
            }
            console.log(round);
        }

        console.log(capacityNeverExceeded);
        console.log("Repeat: ", repeat);

        if (strategy[populationIndex] == 2) {
            if (!capacityNeverExceeded || lowContactLevel[populationIndex] > 100) {
                lowContactLevel[populationIndex] -= 1;
                highContactLevel[populationIndex] -= 1;
                lastRound = true;
                initializeValues();
            }
            else {
                lowContactLevel[populationIndex] += 1;
                highContactLevel[populationIndex] += 1;
                lastRound = false;
                capacityNeverExceeded = true;
                round = numberSimulations[populationIndex];
            }
        }
    }

    // determine maximum y value for chart 1
    var maxValueChart = [];
    maxValueChart[0] = 0;
    for (i = 0; i < durationSimulation[populationIndex] + zeroPoint; i++) {
        if (hospitalCapacity[populationIndex][i] > maxValueChart[0]) maxValueChart[0] = hospitalCapacity[populationIndex][i];
        if (totalMaxNeedingIntensiveCare[populationIndex][i] > maxValueChart[0]) maxValueChart[0] = totalMaxNeedingIntensiveCare[populationIndex][i];
    }

    // determine maximum y value for chart 4
    maxValueChart[1] = 0;
    for (i = 0; i < durationSimulation[populationIndex] + zeroPoint; i++) {
        if (totalMaxDeceased[populationIndex][i] > maxValueChart[1]) maxValueChart[1] = totalMaxDeceased[populationIndex][i];
    }

    for (i = 0; i < 2; i++) {
        if (Math.log10(maxValueChart[i]) > 0) maxValueChart[i] = Math.ceil(maxValueChart[i] / Math.pow(10, Math.floor(Math.log10(maxValueChart[i])))) * Math.pow(10, (Math.floor(Math.log10(maxValueChart[i]))));
    }




    // extract data for overview (hospital capacity)
    var viewData1 = [];
    for (i = 0; i < durationSimulation[0] + zeroPoint + 1; i++) {
        viewData1[i] = [];
    }

    for (i = 0; i < 5 * numberPopulations + 1; i++)
        viewData1[0][i] = "";

    viewData1[0][0] = "time";
    viewData1[0][1] = "";
    viewData1[0][2] = legendRequiredIntensiveCareUnits;
    viewData1[0][3] = "";
    viewData1[0][4] = legendAvailableIntensiveCareUnits;
    viewData1[0][5] = "hospitalized";

    for (i = 1; i < durationSimulation[populationIndex] + zeroPoint + 1; i++) {
        viewData1[i][0] = i - zeroPoint;
        for (j = 0; j < numberPopulations; j++) {
            viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 1] = totalMaxNeedingIntensiveCare[j][i];
            viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 2] = totalNeedingIntensiveCare[j][i];
            viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 3] = totalMinNeedingIntensiveCare[j][i];
            viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 4] = hospitalCapacity[j][i];
            viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 5] = totalHospitalized[j][i];
        }
    }

    // produce stacked diagram
    for (i = 1; i < durationSimulation[populationIndex] + zeroPoint + 1; i++) {
        for (j = populationIndex; j < numberPopulations + populationIndex; j++) {
            if (i > 0) {
                viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 5] = viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 5] - viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 4];
                viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 4] = viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 4] - viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 3];
                viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 3] = viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 3] - viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 2];
                viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 2] = viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 2] - viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 1];


                if (j > populationIndex) {
                    viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 1] -= viewData1[i][5 * ((j - 1 + populationIndex) % numberPopulations) + 5];
                    viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 1] -= viewData1[i][5 * ((j - 1 + populationIndex) % numberPopulations) + 4];
                    viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 1] -= viewData1[i][5 * ((j - 1 + populationIndex) % numberPopulations) + 3];
                    viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 1] -= viewData1[i][5 * ((j - 1 + populationIndex) % numberPopulations) + 2];
                    viewData1[i][5 * ((j + populationIndex) % numberPopulations) + 1] -= viewData1[i][5 * ((j - 1 + populationIndex) % numberPopulations) + 1];
                }

            }
        }
    }

    // extract data for overview (infection and immunity)
    var viewData2 = [];
    for (i = 0; i < durationSimulation[0] + zeroPoint + 1; i++) {
        viewData2[i] = [];
    }

    for (i = 0; i < 3 * numberPopulations + 1; i++)
        viewData2[0][i] = "";

    viewData2[0][0] = "time";
    viewData2[0][1] = legendS;
    viewData2[0][2] = legendI;
    viewData2[0][3] = legendR;


    for (i = 1; i < durationSimulation[populationIndex] + zeroPoint + 1; i++) {
        viewData2[i][0] = i - zeroPoint;
        for (j = 0; j < numberPopulations; j++) {
            viewData2[i][3 * ((j + populationIndex) % numberPopulations) + 1] = totalSusceptible[j % numberPopulations][i];
            viewData2[i][3 * ((j + populationIndex) % numberPopulations) + 2] = totalInfectious[j % numberPopulations][i];
            viewData2[i][3 * ((j + populationIndex) % numberPopulations) + 3] = totalRemoved[j % numberPopulations][i];
        }
    }

    var viewData3 = [];
    for (i = 0; i < durationSimulation[populationIndex] + zeroPoint + 1; i++) {
        viewData3[i] = [];
    }

    for (i = 0; i < 3 * numberPopulations + 1; i++)
        viewData3[0][i] = "";

    viewData3[0][0] = "time";
    viewData3[0][1] = legendDesiredContactReduction;
    viewData3[0][2] = legendEffectiveContactReduction;
    viewData3[0][3] = legendAverageEffectiveContactReduction;


    for (i = 1; i < durationSimulation[populationIndex] + zeroPoint + 1; i++) {
        viewData3[i][0] = i - zeroPoint;
        for (j = 0; j < numberPopulations; j++) {
            viewData3[i][3 * ((j + populationIndex) % numberPopulations) + 1] = desiredContactReduction[j % numberPopulations][i];
            viewData3[i][3 * ((j + populationIndex) % numberPopulations) + 2] = actualContactReduction[j % numberPopulations][i];
            viewData3[i][3 * ((j + populationIndex) % numberPopulations) + 3] = averageContactReduction[j % numberPopulations][i];

            if (i > 0) {
                viewData3[i][3 * ((j + populationIndex) % numberPopulations) + 1] = viewData3[i][3 * ((j + populationIndex) % numberPopulations) + 1] / 100;
                viewData3[i][3 * ((j + populationIndex) % numberPopulations) + 2] = viewData3[i][3 * ((j + populationIndex) % numberPopulations) + 2] / 100;
                viewData3[i][3 * ((j + populationIndex) % numberPopulations) + 3] = viewData3[i][3 * ((j + populationIndex) % numberPopulations) + 3] / 100;
            }
        }
    }

    var viewData4 = [];
    for (i = 0; i < durationSimulation[populationIndex] + zeroPoint + 1; i++) {
        viewData4[i] = [];
    }

    for (i = 0; i < 3 * numberPopulations + 1; i++)
        viewData4[0][i] = "";

    viewData4[0][0] = "time";
    viewData4[0][1] = "";
    viewData4[0][2] = legendDeceased;
    viewData4[0][3] = "";

    for (i = 0; i < durationSimulation[populationIndex] + zeroPoint + 0; i++) {
        viewData4[i + 1][0] = i - zeroPoint;
        for (j = 0; j < numberPopulations; j++) {
            viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 3] = totalMaxDeceased[j][i];
            viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 2] = totalDeceased[j][i];
            viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 1] = totalMinDeceased[j][i];
        }
        for (j = populationIndex; j < numberPopulations + populationIndex; j++) {
            // produce stacked diagram
            if (i >= 0) {
                viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 3] = viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 3] - viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 2];
                viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 2] = viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 2] - viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 1];
                viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 1] = viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 1];

                // HERE
                if (j > populationIndex) {
                    viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 1] -= viewData4[i + 1][3 * ((j - 1 + populationIndex) % numberPopulations) + 3];
                    viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 1] -= viewData4[i + 1][3 * ((j - 1 + populationIndex) % numberPopulations) + 2];
                    viewData4[i + 1][3 * ((j + populationIndex) % numberPopulations) + 1] -= viewData4[i + 1][3 * ((j - 1 + populationIndex) % numberPopulations) + 1];
                }

            }
        }
    }

    // graphical representation

    google.charts.load('current', { 'packages': ['corechart'] });
    google.charts.setOnLoadCallback(drawChart);

    document.body.style.cursor = 'default';
    var panels = document.getElementsByTagName("*");
    for (i = 0; i < panels.length; i++) {
        panels[i].style.cursor = 'default';
    }



    // functions

    function calculateContactMatrix(index, optionChoice) {

        // define variables
        var i, j, l, m;                         // counters
        var normingFactors = [];                // factors used to norm the contactMatrix such that cross-border interaction is compensated (same amount of contact no matter if the border is closed or not)

        // build 2 dimensional array
        for (i = 0; i < numberPopulations; i++) {
            normingFactors[i] = [];
        }

        for (j = 0; j < numberSubPopulations; j++) {
            for (l = 0; l < numberPopulations; l++) {
                for (m = 0; m < numberSubPopulations; m++) {
                    if (index == l) {
                        contactMatrix[index][j][l][m] = subpopulation[l][m] / population[l] * (population[l] - borderCrossers[index]) / population[l];
                    }
                    else {
                        contactMatrix[index][j][l][m] = borderCrossers[index] / population[l] * subpopulation[l][m] / population[l];
                    }
                    if (optionChoice == 1) {
                        if (j == 1) contactMatrix[index][j][l][m] = 0;
                        if (m == 1) contactMatrix[index][j][l][m] = 0;
                    }
                }
            }
        }
        // determine norming factors for contact matrix
        for (j = 0; j < numberSubPopulations; j++) {
            normingFactors[index][j] = 0;
            for (l = 0; l < numberPopulations; l++) {
                for (m = 0; m < numberSubPopulations; m++) {
                    normingFactors[index][j] += contactMatrix[index][j][l][m];
                }
            }
        }
        // norm contact matrix
        for (j = 0; j < numberSubPopulations; j++) {
            for (l = 0; l < numberPopulations; l++) {
                for (m = 0; m < numberSubPopulations; m++) {
                    if (normingFactors[index][j] == 0) contactMatrix[index][j][l][m] = 0;
                    else contactMatrix[index][j][l][m] /= normingFactors[index][j];
                }
            }
        }
    }


    function drawChart() {

        // visualization of Figure A

        var data1 = google.visualization.arrayToDataTable(viewData1);
        var options1 = {
            title: captionFigureA,
            curveType: 'none',
            fontSize: graphFontSize,
            legend: { position: 'bottom' },
            vAxes: {
                0: {
                    format: 'short',
                    viewWindowMode: 'explicit',
                    viewWindow: {
                        max: maxValueChart[0],
                        min: 0,
                    },
                },
                1: { format: 'short' },
            },
            isStacked: 'true',
            hAxis: {
                title: legendTime,
            },
            seriesType: 'area',
            enableInteractivity: false,
            series: {
                0: { targetAxisIndex: 0, color: '#3366CC', lineWidth: 0, areaOpacity: 0, visibleInLegend: false, interpolateNulls: false },
                1: { targetAxisIndex: 0, color: '#3366CC', lineWidth: 2, areaOpacity: 0.2, interpolateNulls: false },
                2: { targetAxisIndex: 0, color: '#3366CC', lineWidth: 0, areaOpacity: 0.2, visibleInLegend: false, interpolateNulls: false },
                3: { targetAxisIndex: 0, color: '#DC3912', lineWidth: 2, areaOpacity: 0, interpolateNulls: false },
                4: { targetAxisIndex: 0, color: '#FF9900', pointShape: 'square', pointSize: 0, lineWidth: 0, visibleInLegend: false, areaOpacity: 0, interpolateNulls: false },

                5: { targetAxisIndex: 0, color: '#b8c7e6', lineDashStyle: [6, 6], lineWidth: 0, areaOpacity: 0, visibleInLegend: false, interpolateNulls: false },
                6: { targetAxisIndex: 0, color: '#b8c7e6', lineDashStyle: [6, 6], lineWidth: 2, areaOpacity: 0.2, visibleInLegend: false, interpolateNulls: false },
                7: { targetAxisIndex: 0, color: '#b8c7e6', lineDashStyle: [6, 6], lineWidth: 0, areaOpacity: 0.2, visibleInLegend: false, interpolateNulls: false },
                8: { targetAxisIndex: 0, color: '#dba79a', lineDashStyle: [6, 6], lineWidth: 2, areaOpacity: 0, visibleInLegend: false, interpolateNulls: false },
                9: { targetAxisIndex: 0, color: '#ffe1b3', lineDashStyle: [6, 6], lineWidth: 0, areaOpacity: 0, visibleInLegend: false, interpolateNulls: false },
            },
            backgroundColor: 'white',
            backgroundColor: { stroke: 'black', strokeWidth: 2 },
        };

        var chart1 = new google.visualization.ComboChart(document.getElementById('curve_chart_0'));
        chart1.draw(data1, options1);


        // visualization of Figure B

        var data2 = google.visualization.arrayToDataTable(viewData2);

        // add annotations
        data2.insertColumn(4, { type: 'string', role: 'annotationText' });
        data2.insertColumn(4, { type: 'string', role: 'annotation' });

        data2.insertColumn(3, { type: 'string', role: 'annotationText' });
        data2.insertColumn(3, { type: 'string', role: 'annotation' });

        data2.insertColumn(2, { type: 'string', role: 'annotationText' });
        data2.insertColumn(2, { type: 'string', role: 'annotation' });

        data2.addColumn({ type: 'string', role: 'annotation' });
        data2.addColumn({ type: 'string', role: 'annotationText' });


        if (option[populationIndex] == 1) {
            data2.setCell(1, 2, 'Empfängliche (susceptible)', noteRiskGroupSeparated);
            data2.setCell(durationIntervention[populationIndex], 2, 'Empfängliche (susceptible)', noteRiskGroupIntegrated);
        }

        if (herdImmunityAlwaysReached[populationIndex]) data2.setCell(iHerdImmunity[populationIndex], 8, 'Genesene (recovered)', noteHerdImmunityAlwaysReached);
        else data2.setCell(durationSimulation[populationIndex] + zeroPoint - 2, 8, 'Genesene (recovered)', noteHerdImmunityNotAlwaysReached);

        if (startLevelAlwaysReached[populationIndex]) data2.setCell(iStartLevelReached[populationIndex], 5, 'Infizierte (infected)', noteRiskAlwaysReduced);
        else data2.setCell(durationSimulation[populationIndex] + zeroPoint - 2, 5, 'Infizierte (infected)', noteRiskNotAlwaysReduced);

        var options2 = {
            title: captionFigureB,
            curveType: 'none',
            fontSize: graphFontSize,
            legend: { position: 'bottom' },
            hAxis: {
                title: legendTime,
            },
            vAxes: {
                0: { format: 'short', type: "line", targetAxisIndex: 0 },
                1: { format: 'short', type: "line", targetAxisIndex: 1 },
            },
            series: {
                0: { targetAxisIndex: 0, color: '#3366CC', annotations: { style: { type: 'line', }, stem: { length: 50, color: '#3366CC' } }, interpolateNulls: false },
                1: { targetAxisIndex: 0, color: '#DC3912', annotations: { style: { type: 'line', }, stem: { length: 20, color: '#DC3912' } }, interpolateNulls: false },
                2: { targetAxisIndex: 0, color: '#FF9900', annotations: { style: { type: 'line', }, stem: { length: 50, color: '#FF9900' } }, interpolateNulls: false },
                3: { targetAxisIndex: 0, color: '#b8c7e6', lineDashStyle: [6, 6], annotations: { style: { type: 'line', }, stem: { length: 50, color: '#a1b8e6' } }, interpolateNulls: false, visibleInLegend: false },
                4: { targetAxisIndex: 0, color: '#dba79a', lineDashStyle: [6, 6], annotations: { style: { type: 'line', }, stem: { length: 20, color: '#a1b8e6' } }, interpolateNulls: false, visibleInLegend: false },
                5: { targetAxisIndex: 0, color: '#ffe1b3', lineDashStyle: [6, 6], annotations: { style: { type: 'line', }, stem: { length: 50, color: '#a1b8e6' } }, interpolateNulls: false, visibleInLegend: false },

            },
            backgroundColor: { stroke: 'black', strokeWidth: 2 },
        };

        var chart2 = new google.visualization.ComboChart(document.getElementById('curve_chart_1'));
        chart2.draw(data2, options2);

        // visualization of Figure C

        var data3 = google.visualization.arrayToDataTable(viewData3);
        var options3 = {
            title: captionFigureC,
            curveType: 'none',
            fontSize: graphFontSize,
            legend: { position: 'bottom' },
            vAxes: {
                0: { format: '#.###%' }
            },
            hAxis: {
                title: legendTime,
            },
            series: {
                1: { targetAxisIndex: 0, color: '#3366CC', lineWidth: 2, areaOpacity: 0 },
                2: { targetAxisIndex: 0, color: '#DC3912', lineWidth: 2, areaOpacity: 0 },
                0: { targetAxisIndex: 0, color: '#FF9900', lineWidth: 2, areaOpacity: 0 },
                4: { targetAxisIndex: 0, color: '#b8c7e6', lineDashStyle: [6, 6], annotations: { style: { type: 'line', }, stem: { length: 50, color: '#a1b8e6' } }, interpolateNulls: false, visibleInLegend: false },
                5: { targetAxisIndex: 0, color: '#dba79a', lineDashStyle: [6, 6], annotations: { style: { type: 'line', }, stem: { length: 20, color: '#a1b8e6' } }, interpolateNulls: false, visibleInLegend: false },
                3: { targetAxisIndex: 0, color: '#ffe1b3', lineDashStyle: [6, 6], annotations: { style: { type: 'line', }, stem: { length: 50, color: '#a1b8e6' } }, interpolateNulls: false, visibleInLegend: false },

            },
            vAxis: {
                viewWindowMode: 'explicit',
                viewWindow: {
                    max: 1,
                    min: 0
                }
            },
            backgroundColor: { stroke: 'black', strokeWidth: 2 },
        };
        var chart3 = new google.visualization.ComboChart(document.getElementById('curve_chart_2'));
        chart3.draw(data3, options3);

        // visualization of Figure D
        var data4 = google.visualization.arrayToDataTable(viewData4);
        var options4 = {
            title: captionFigureD,
            curveType: 'none',
            fontSize: graphFontSize,
            legend: { position: 'bottom' },
            vAxes: {
                0: {
                    format: 'short',
                    viewWindowMode: 'explicit',
                    viewWindow: {
                        max: maxValueChart[1],
                        min: 0,
                    },
                },
                1: { format: 'short' }
            },
            isStacked: 'true',
            hAxis: {
                title: legendTime,
            },
            seriesType: 'area',
            enableInteractivity: false,
            series: {
                0: { targetAxisIndex: 0, color: '#3366CC', lineWidth: 0, areaOpacity: 0, visibleInLegend: false },
                1: { targetAxisIndex: 0, color: '#3366CC', lineWidth: 2, areaOpacity: 0.2 },
                2: { targetAxisIndex: 0, color: '#3366CC', lineWidth: 0, areaOpacity: 0.2, visibleInLegend: false },
                3: { targetAxisIndex: 0, color: '#b8c7e6', lineDashStyle: [6, 6], lineWidth: 0, areaOpacity: 0, visibleInLegend: false },
                4: { targetAxisIndex: 0, color: '#b8c7e6', lineDashStyle: [6, 6], lineWidth: 2, areaOpacity: 0.2, visibleInLegend: false },
                5: { targetAxisIndex: 0, color: '#b8c7e6', lineDashStyle: [6, 6], lineWidth: 0, areaOpacity: 0.2, visibleInLegend: false },
            },
            backgroundColor: { stroke: 'black', strokeWidth: 2 },
        };
        var chart4 = new google.visualization.ComboChart(document.getElementById('curve_chart_3'));
        chart4.draw(data4, options4);

    }

    // export data is requested
    if (exportData) {
        var writeData = [];
        writeData[0] = "\ufeff\n";

        // write headlines
        writeData[1] = "totalMinNeedingIntensiveCare[0]";
        writeData[2] = "\t";
        writeData[3] = "totalNeedingIntensiveCare[0]";
        writeData[4] = "\t";
        writeData[5] = "totalMaxNeedingIntensiveCare[0]";
        writeData[6] = "\t";
        writeData[7] = "hospitalCapacity[0]";
        writeData[8] = "\t";
        writeData[9] = "totalSusceptible[0]";
        writeData[10] = "\t";
        writeData[11] = "totalInfectious[0]";
        writeData[12] = "\t";
        writeData[13] = "totalRemoved[0]";
        writeData[14] = "\t";
        writeData[15] = "desiredContactReduction[0]";
        writeData[16] = "\t";
        writeData[17] = "actualContactReduction[0]";
        writeData[18] = "\t";
        writeData[19] = "averageContactReduction[0]";
        writeData[20] = "\t";
        writeData[21] = "totalMinDeceased[0]";
        writeData[22] = "\t";
        writeData[23] = "totalDeceased[0]";
        writeData[24] = "\t";
        writeData[25] = "totalMaxDeceased[0]";
        writeData[26] = "\t";

        if (numberPopulations > 1) {
            writeData[27] = "totalMinNeedingIntensiveCare[1]";
            writeData[28] = "\t";
            writeData[29] = "totalNeedingIntensiveCare[1]";
            writeData[30] = "\t";
            writeData[31] = "totalMaxNeedingIntensiveCare[1]";
            writeData[32] = "\t";
            writeData[33] = "hospitalCapacity[1]";
            writeData[34] = "\t";
            writeData[35] = "totalSusceptible[1]";
            writeData[36] = "\t";
            writeData[37] = "totalInfectious[1]";
            writeData[38] = "\t";
            writeData[39] = "totalRemoved[1]";
            writeData[40] = "\t";
            writeData[41] = "desiredContactReduction[1]";
            writeData[42] = "\t";
            writeData[43] = "actualContactReduction[1]";
            writeData[44] = "\t";
            writeData[45] = "averageContactReduction[1]";
            writeData[46] = "\t";
            writeData[47] = "totalMinDeceased[1]";
            writeData[48] = "\t";
            writeData[49] = "totalDeceased[1]";
            writeData[50] = "\t";
            writeData[51] = "totalMaxDeceased[1]";
            writeData[52] = "\n";
        }

        // write data
        for (i = 0; i < durationSimulation[0] + zeroPoint; i++) {
            writeData[26 * numberPopulations * (i + 1) + 1] = totalMinNeedingIntensiveCare[0][i];
            writeData[26 * numberPopulations * (i + 1) + 2] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 3] = totalNeedingIntensiveCare[0][i];
            writeData[26 * numberPopulations * (i + 1) + 4] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 5] = totalMaxNeedingIntensiveCare[0][i];
            writeData[26 * numberPopulations * (i + 1) + 6] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 7] = hospitalCapacity[0][i];
            writeData[26 * numberPopulations * (i + 1) + 8] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 9] = totalSusceptible[0][i];
            writeData[26 * numberPopulations * (i + 1) + 10] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 11] = totalInfectious[0][i];
            writeData[26 * numberPopulations * (i + 1) + 12] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 13] = totalRemoved[0][i];
            writeData[26 * numberPopulations * (i + 1) + 14] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 15] = desiredContactReduction[0][i];
            writeData[26 * numberPopulations * (i + 1) + 16] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 17] = actualContactReduction[0][i];
            writeData[26 * numberPopulations * (i + 1) + 18] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 19] = averageContactReduction[0][i];
            writeData[26 * numberPopulations * (i + 1) + 20] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 21] = totalMinDeceased[0][i];
            writeData[26 * numberPopulations * (i + 1) + 22] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 23] = totalDeceased[0][i];
            writeData[26 * numberPopulations * (i + 1) + 24] = "\t";
            writeData[26 * numberPopulations * (i + 1) + 25] = totalMaxDeceased[0][i];
            writeData[26 * numberPopulations * (i + 1) + 26] = "\t";

            if (numberPopulations > 1) {
                writeData[26 * numberPopulations * (i + 1) + 27] = totalMinNeedingIntensiveCare[1][i];
                writeData[26 * numberPopulations * (i + 1) + 28] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 29] = totalNeedingIntensiveCare[1][i];
                writeData[26 * numberPopulations * (i + 1) + 30] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 31] = totalMaxNeedingIntensiveCare[1][i];
                writeData[26 * numberPopulations * (i + 1) + 32] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 33] = hospitalCapacity[1][i];
                writeData[26 * numberPopulations * (i + 1) + 34] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 35] = totalSusceptible[1][i];
                writeData[26 * numberPopulations * (i + 1) + 36] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 37] = totalInfectious[1][i];
                writeData[26 * numberPopulations * (i + 1) + 38] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 39] = totalRemoved[1][i];
                writeData[26 * numberPopulations * (i + 1) + 40] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 41] = desiredContactReduction[1][i];
                writeData[26 * numberPopulations * (i + 1) + 42] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 43] = actualContactReduction[1][i];
                writeData[26 * numberPopulations * (i + 1) + 44] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 45] = averageContactReduction[1][i];
                writeData[26 * numberPopulations * (i + 1) + 46] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 47] = totalMinDeceased[1][i];
                writeData[26 * numberPopulations * (i + 1) + 48] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 49] = totalDeceased[1][i];
                writeData[26 * numberPopulations * (i + 1) + 50] = "\t";
                writeData[26 * numberPopulations * (i + 1) + 51] = totalMaxDeceased[1][i];
                writeData[26 * numberPopulations * (i + 1) + 52] = "\n";
            }
        }
        csv = new Blob(writeData, { type: "text/csv;charset=utf-8", endings: 'native' });

        var link = document.createElement('a');
        link.style.display = "";
        link.setAttribute('href', window.URL.createObjectURL(csv));
        link.setAttribute('download', 'export.csv');
        document.body.append(link);
        link.click();
        document.body.removeChild(link);
    }
}