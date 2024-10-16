#include "LKH.h"
#include "Genetic.h"

/*
 * This file contains the main function of the program.
 */

int LKHmain(double startTime)
{
    GainType Cost, OldOptimum;
    double Time, LastTime;

    StartTime = startTime;
    LastTime = GetTime();

    AllocateStructures();
    CreateCandidateSet();
    InitializeStatistics();

    if (Norm != 0)
        BestCost = PLUS_INFINITY;
    else {
        /* The ascent has solved the problem! */
        Optimum = BestCost = (GainType) LowerBound;
        UpdateStatistics(Optimum, GetTime() - LastTime);
        RecordBetterTour();
        RecordBestTour();
        Runs = 0;
    }

    /* Find a specified number (Runs) of local optima */
    for (Run = 1; Run <= Runs; Run++) {
        LastTime = GetTime();
        // if (Run!=1 && LastTime - StartTime >= TotalTimeLimit) {
        //     if (TraceLevel >= 1)
        //         printff("*** Time limit exceeded in LKHmain ***\n");
        //     Run--;
        //     break;
        // }
        Cost = FindTour();      /* using the Lin-Kernighan heuristic */
        
        if (Cost < BestCost) {
            BestCost = Cost;
            RecordBetterTour();
            RecordBestTour();
        }
        OldOptimum = Optimum;
        if (Cost < Optimum) {
            if (FirstNode->InputSuc) {
                Node *N = FirstNode;
                while ((N = N->InputSuc = N->Suc) != FirstNode);
            }
            Optimum = Cost;
            printff("*** New optimum = " GainFormat " ***\n", Optimum);
        }
        Time = fabs(GetTime() - LastTime);
        UpdateStatistics(Cost, Time);
        if (TraceLevel >= 1 && Cost != PLUS_INFINITY) {
            printff("Run %d: Cost = " GainFormat, Run, Cost);
            if (Optimum != MINUS_INFINITY && Optimum != 0)
                printff(", Gap = %0.4f%%",
                        100.0 * (Cost - Optimum) / Optimum);
            printff(", Time = %0.2f sec. %s\n\n", Time,
                    Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
        }
        SRandom(++Seed);
    }
    PrintStatistics();
    return EXIT_SUCCESS;
}
