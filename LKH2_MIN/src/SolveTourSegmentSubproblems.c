#include "LKH.h"

/*
 * The SolveTourSegmentSubproblems function breaks the tour up into 
 * segments and attempts to find an improved subtour for each segment. 
 * If an improvement is found, the new tour is written to TourFile.
 * 
 * The tour is divided into segments of equal size. Each segment 
 * constitutes a subproblem to be solved by the LinKernighan function. 
 * 
 * The original tour is given by the SubproblemSuc references of the nodes. 
 * The size of each segment is SubproblemSize.
 *
 * The function runs in two rounds. In the first round the first segment 
 * starts at FirstNode. In the second round the middle of this segment is 
 * used as a starting point for the first segment. 
 */

void SolveTourSegmentSubproblems()
{
    Node *FirstNodeSaved, *N;
    int CurrentSubproblem, Subproblems, Round, i;
    GainType GlobalBestCost, OldGlobalBestCost;
    double EntryTime = GetTime();

    AllocateStructures();
    Subproblems = (int) ceil((double) Dimension / SubproblemSize);
    ReadPenalties();
    FirstNode = &NodeSet[Random() % Dimension + 1];

    /* Compute upper bound for the original problem */
    GlobalBestCost = 0;
    N = FirstNodeSaved = FirstNode;
    do {
        if (!Fixed(N, N->SubproblemSuc))
            GlobalBestCost += Distance(N, N->SubproblemSuc);
        N->Subproblem = 0;
    }
    while ((N = N->SubproblemSuc) != FirstNode);

    // for (Round = 1; Round <= 2; Round++) {
    for (Round = 1; Round <= 1; Round++) { // 只一轮
        if (Round == 2 && Subproblems == 1)
            break;
        if (TraceLevel >= 1) {
            if (Round == 2 || TraceLevel >= 2)
                printff("\n");
            printff
                ("*** Tour segment partitioning *** [Round %d of %d, Cost = "
                //  GainFormat "]\n", Round, Subproblems > 1 ? 2 : 1,
                 GainFormat "]\n", Round, 1,
                 GlobalBestCost);
        }
        FirstNode = FirstNodeSaved;
        if (Round == 2)
            for (i = SubproblemSize / 2; i > 0; i--)
                FirstNode = FirstNode->SubproblemSuc;

        /* 求解每个子问题 */
        for (CurrentSubproblem = 1; CurrentSubproblem <= Subproblems; CurrentSubproblem++) {
            // 计算当前子任务的总时间限制
            SubProblemStartTime = GetTime();
            SubProblemTotalTimeLimit = (TotalTimeLimit - (SubProblemStartTime - StartTime)) / (Subproblems-CurrentSubproblem+1);
            for (i = 0, N = FirstNode;
                 i < SubproblemSize ||
                 (FirstNode->Id <= DimensionSaved) != (N->Id <=
                                                       DimensionSaved);
                 i++, N = N->SubproblemSuc) {
                N->Subproblem =
                    (Round - 1) * Subproblems + CurrentSubproblem;
                N->FixedTo1Saved = N->FixedTo2Saved = 0;
                N->SubBestPred = N->SubBestSuc = 0;
            }
            OldGlobalBestCost = GlobalBestCost;
            SolveSubproblem((Round - 1) * Subproblems + CurrentSubproblem,
                            Subproblems, &GlobalBestCost);
            // if (SubproblemsCompressed
            //     && GlobalBestCost == OldGlobalBestCost)
            //     SolveCompressedSubproblem((Round - 1) * Subproblems +
            //                               CurrentSubproblem, Subproblems,
            //                               &GlobalBestCost);
            FirstNode = N;
        }
    }
    printff("\nCost = " GainFormat, GlobalBestCost);
    if (Optimum != MINUS_INFINITY && Optimum != 0)
        printff(", Gap = %0.4f%%",
                100.0 * (GlobalBestCost - Optimum) / Optimum);
    printff(", Time = %0.2f sec. %s\n", fabs(GetTime() - EntryTime),
            GlobalBestCost < Optimum ? "<" : GlobalBestCost ==
            Optimum ? "=" : "");
    // if (SubproblemBorders && Subproblems > 1)
    //     SolveSubproblemBorderProblems(Subproblems, &GlobalBestCost);
}
