#include "LKH.h"

/*
 * After the candidate set has been created the FindTour function is called
 * a predetermined number of times (Runs).
 *
 * FindTour performs a number of trials, where in each trial it attempts
 * to improve a chosen initial tour using the modified Lin-Kernighan edge
 * exchange heuristics.
 *
 * Each time a better tour is found, the tour is recorded, and the candidates
 * are reorderded by the AdjustCandidateSet function. Precedence is given to
 * edges that are common to two currently best tours. The candidate set is
 * extended with those tour edges that are not present in the current set.
 * The original candidate set is re-established at exit from FindTour.
 */

static void SwapCandidateSets();
static GainType OrdinalTourCost;

GainType FindTour()
{
    GainType Cost;
    Node *t;
    int i;
    double EntryTime = GetTime();

    t = FirstNode;
    do
        t->OldPred = t->OldSuc = t->NextBestSuc = t->BestSuc = 0;
    while ((t = t->Suc) != FirstNode);
    if (Run == 1 && Dimension == DimensionSaved) {
        OrdinalTourCost = 0;
        for (i = 1; i < Dimension; i++)
            OrdinalTourCost += C(&NodeSet[i], &NodeSet[i + 1])
                - NodeSet[i].Pi - NodeSet[i + 1].Pi;
        OrdinalTourCost += C(&NodeSet[Dimension], &NodeSet[1])
            - NodeSet[Dimension].Pi - NodeSet[1].Pi;
        OrdinalTourCost /= Precision;
    }
    BetterCost = PLUS_INFINITY;
    if (MaxTrials > 0)
        HashInitialize(HTable);
    else {
        Trial = 1;
        ChooseInitialTour();
    }

    double recordTime = GetTime();
    GainType recordCost = BetterCost;
    // int recordTrial = 1;
    // int TrialSpan = 100;
    // if(DimensionSaved <= 1002)
    //     TrialSpan = 400;
    // else if(DimensionSaved <= 2002)
    //     TrialSpan = 58;
    // else if(DimensionSaved <= 5002)
    //     TrialSpan = 12;
    // else
    //     TrialSpan = 4;

    for (Trial = 1; Trial <= MaxTrials; Trial++)
    {  
        // 总时间限制
        if (GetTime() - StartTime >= TotalTimeLimit) {
            if (TraceLevel >= 1)
                printff("*** Time limit exceeded in FindTour ***\n");
            Trial--;
            break;
        }

        // 在一定时间跨度内统计改进幅度
        if(GetTime() - recordTime >= TimeSpan){
            // 20s之前
            if(GetTime() - StartTime <= 20 && recordCost - BetterCost < TimeSpan*ScheduleScoreInSecond){
                if (TraceLevel >= 1)
                    printff("*** Within 20s the extent of improvement("GainFormat"(<%.1f)) is too small in %.1fs ***\n",(recordCost - BetterCost),TimeSpan*ScheduleScoreInSecond, TimeSpan);
                Trial--;
                break;
            }
            // 20s之后
            else if(GetTime() - StartTime > 20 && recordCost - BetterCost < TimeSpan*PenaltyScoreInSecond){
                if (TraceLevel >= 1)
                    printff("*** The extent of improvement("GainFormat"(<%.1f)) is too small in %.1fs ***\n",(recordCost - BetterCost),TimeSpan*PenaltyScoreInSecond, TimeSpan);
                Trial--;
                break;
            }
            else{
                recordTime = GetTime();
                recordCost = BetterCost;
            }
        }

        // 在一定Trial跨度内统计改进幅度
        // if(Trial - recordTrial >= TimeSpan * TrialSpan){
        //     if(recordCost - BetterCost < TimeSpan*ScheduleScoreInSecond){
        //         if (TraceLevel >= 1)
        //             printff("*** The extent of improvement("GainFormat") is too small in %d Trials ***\n",(recordCost - BetterCost), (int)TimeSpan * TrialSpan);
        //         Trial--;
        //         break;
        //     }
        //     else{
        //         recordTrial = Trial;
        //         recordCost = BetterCost;
        //     }
        // }

        /* Choose FirstNode at random */
        if (Dimension == DimensionSaved)
            FirstNode = &NodeSet[1 + Random() % Dimension];
        else
            for (i = Random() % Dimension; i > 0; i--)
                FirstNode = FirstNode->Suc;
        ChooseInitialTour();
        Cost = LinKernighan();
        // if (Trial == 1 || (GetTime() - EntryTime < TimeLimit &&
        //                    GetTime() - StartTime < TotalTimeLimit)) {
        if (FirstNode->BestSuc) {
            /* Merge tour with current best tour */
            t = FirstNode;
            while ((t = t->Next = t->BestSuc) != FirstNode);
            Cost = MergeWithTour();
        }
        if (Dimension == DimensionSaved && Cost >= OrdinalTourCost &&
            BetterCost > OrdinalTourCost) {
            /* Merge tour with ordinal tour */
            for (i = 1; i < Dimension; i++)
                NodeSet[i].Next = &NodeSet[i + 1];
            NodeSet[Dimension].Next = &NodeSet[1];
            Cost = MergeWithTour();
        }
        // }
        if (Cost < BetterCost) {
            if (TraceLevel >= 1) {
                printff("* %d: Cost = " GainFormat, Trial, Cost);
                if (Optimum != MINUS_INFINITY && Optimum != 0)
                    printff(", Gap = %0.4f%%",
                            100.0 * (Cost - Optimum) / Optimum);
                printff(", Time = %0.2f sec. %s\n",
                        fabs(GetTime() - EntryTime),
                        Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
            }
            BetterCost = Cost;
            RecordBetterTour();

            if (StopAtOptimum && BetterCost == Optimum)
                break;
            AdjustCandidateSet();
            HashInitialize(HTable);
            HashInsert(HTable, Hash, Cost);
        } else if (TraceLevel >= 2)
            printff("  %d: Cost = " GainFormat ", Time = %0.2f sec.\n",
                    Trial, Cost, fabs(GetTime() - EntryTime));
        /* Record backbones if wanted */
        if (Trial <= BackboneTrials && BackboneTrials < MaxTrials) {
            SwapCandidateSets();
            AdjustCandidateSet();
            if (Trial == BackboneTrials) {
                if (TraceLevel >= 1) {
                    printff("# %d: Backbone candidates ->\n", Trial);
                    CandidateReport();
                }
            } else
                SwapCandidateSets();
        }
    }
    if (BackboneTrials > 0 && BackboneTrials < MaxTrials) {
        if (Trial > BackboneTrials ||
            (Trial == BackboneTrials &&
             (!StopAtOptimum || BetterCost != Optimum)))
            SwapCandidateSets();
        t = FirstNode;
        do {
            free(t->BackboneCandidateSet);
            t->BackboneCandidateSet = 0;
        } while ((t = t->Suc) != FirstNode);
    }
    t = FirstNode;
    if (Norm == 0 || MaxTrials == 0 || !t->BestSuc) {
        do
            t = t->BestSuc = t->Suc;
        while (t != FirstNode);
    }
    Hash = 0;
    do {
        (t->Suc = t->BestSuc)->Pred = t;
        Hash ^= Rand[t->Id] * Rand[t->Suc->Id];
    } while ((t = t->BestSuc) != FirstNode);
    if (Trial > MaxTrials)
        Trial = MaxTrials;
    ResetCandidateSet();
    return BetterCost;
}

/*
 * The SwapCandidateSets function swaps the normal and backbone candidate sets.
 */

static void SwapCandidateSets()
{
    Node *t = FirstNode;
    do {
        Candidate *Temp = t->CandidateSet;
        t->CandidateSet = t->BackboneCandidateSet;
        t->BackboneCandidateSet = Temp;
    } while ((t = t->Suc) != FirstNode);
}
