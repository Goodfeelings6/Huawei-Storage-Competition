#include "LKH.h"
#include "Genetic.h"
#include "Heap.h"
#include "LKHInterface.h"

/*
 * LKH内核 接口，欲调用LHK，按以下三个步骤走：
 * (1)设定 LKH内核 的相关参数值：函数 ReadParameters()。一般情况下无需手动设置，采用默认值即可;
 * (2)设定问题类型、规模等，并进行样例数据的输入：函数 ReadProblem()。此处是ATSP问题，邻接矩阵形式输入;
 * (3)调用 LKH内核 求解：LKHmain().
 *
 * 求解完成后：
 * (1)提取求解结果以返回：函数 OutputTourResult() ;
 * (2)重置 LKH内核 状态：函数 ReSetLKH()。这一步确保当前调用不会对后续调用产生影响.
 */

static void CheckSpecificationPart(void);
static char *Copy(char *S);
static void CreateNodes(void);
static int FixEdge(Node * Na, Node * Nb);
static void Read_EDGE_WEIGHT_SECTION(int **adjMat);
static int TwoDWeightType(void);
static int ThreeDWeightType(void);

void ReadParameters()
{
/*
 * The ReadParameters function reads the name of a parameter file from
 * standard input and reads the problem parameters from this file.
 *
 * All entries of the parameter file are of the form <keyword >= <value>
 * (or <keyword><whitespace><value>), where <keyword> denotes an alphanumeric
 * keyword and <value> denotes alphanumeric or numeric data. Keywords are not
 * case sensitive.
 *
 * The order of specifications in the file is arbitrary. The following
 * specification is mandatory.
 *
 * PROBLEM_FILE = <string>
 * Specifies the name of the problem file.
 *
 * Additional control information may be supplied in the following format:
 *
 * ASCENT_CANDIDATES = <integer>
 * The number of candidate edges to be associated with each node during the
 * ascent. The candidate set is complemented such that every candidate edge
 * is associated with both its two end nodes.
 * Default: 50
 *
 * BACKBONE_TRIALS = <integer>
 * The number of backbone trials in each run.
 * Default: 0
 *
 * BACKTRACKING = { YES | NO }
 * Specifies whether a backtracking k-opt move is to be used as the first
 * move in a sequence of moves (where k = MOVE_TYPE).
 * Default: NO
 *
 * BWTSP = <integer> <integer> [ <integer> ]
 * Specifies the three parameters (B, Q, L) to a BWTSP instance.
 *   B: Number of black nodes.
 *   Q: Maximum number of white nodes on "black-to-black" paths.
 *   L: Maximum length of any "black-to-black" path.
 * Default: 0 0
 *
 * CANDIDATE_FILE = <string>
 * Specifies the name of a file to which the candidate sets are to be written.
 * If, however, the file already exists, the candidate edges are read from the
 * file.
 * The first line of the file contains the number of nodes.
 * Each of the following lines contains a node number, the number of
 * the dad of the node in the minimum spanning tree (0, if the node has no dad),
 * the number of candidate edges emanating from the node, followed by the
 * candidate edges.
 * For each candidate edge its end node number and alpha-value are given.
 * It is possible to give more than one CANDIDATE_FILE specification. In this
 * case the given files are read and the union of their candidate edges is
 * used as candidate sets.
 *
 * CANDIDATE_SET_TYPE = { ALPHA | DELAUNAY [ PURE ] | NEAREST-NEIGHBOR |
 *                        POPMUSIC | QUADRANT }
 * Specifies the candidate set type.
 * ALPHA is LKH's default type. ALPHA and POPMUSIC are applicable in general.
 * The other types can only be used for instances given by coordinates.
 * The optional suffix PURE for the DELAUNAY type specifies that only
 * edges of the Delaunay graph are used as candidates.
 * Default: ALPHA
 *
 * COMMENT <string>
 * A comment.
 *
 * # <string>
 * A comment.
 *
 * DEPOT = <integer>
 * Specifies the depot node.
 * Default: 1
 *
 * DISTANCE : <real>
 * The maximum length allowed for each route in a CVRP.
 * Default: DBL_MAX
 *
 * EDGE_FILE = <string>
 * Specifies the name of a file of candidate edges in Concorde format.
 * The first line of the file contains the number of nodes followed by
 * the number of edges.
 * Each of the following lines contains the number of the two end nodes of
 * an edge and its cost. The cost may be omitted. OBS: nodes are numbered
 * from zero.
 * It is possible to give more than one EDGE_FILE specification. In this
 * case the given files are read and the union of their candidate edges is
 * used as candidate sets.
 *
 * EOF
 * Terminates the input data. The entry is optional.
 *
 * EXCESS = <real>
 * The maximum alpha-value allowed for any candidate edge is set to
 * EXCESS times the absolute value of the lower bound of a solution
 * tour (determined by the ascent).
 * Default: 1.0/DIMENSION
 *
 * EXTERNAL_SALESMEN = <integer>
 * Specifies the number of external salesmen for an OCMTSP problem.
 * Default: 0
 *
 * EXTRA_CANDIDATES = <integer> [ SYMMETRIC ]
 * Number of extra candidate edges to be added to the candidate set
 * of each node. Their candidate set type may be specified after the
 * keyword EXTRA_CANDIDATE_SET_TYPE.
 * The integer may be followed by the keyword SYMMETRIC, signifying
 * that these extra candidate edges is to be complemented such
 * that each of them is associated with both its two end nodes.
 * Default: 0
 *
 * EXTRA_CANDIDATE_SET_TYPE = { NEAREST-NEIGHBOR | POPMISIC | QUADRANT }
 * The candidate set type of extra candidate edges.
 * Default: QUADRANT
 *
 * GAIN23 = { YES | NO }
 * Specifies whether the Gain23 function is used.
 * Default: YES
 *
 * GAIN_CRITERION = { YES | NO }
 * Specifies whether Lin and Kernighan's gain criterion is used.
 * Default: YES
 *
 * INITIAL_PERIOD = <integer>
 * The length of the first period in the ascent.
 * Default: DIMENSION/2 (but at least 100)
 *
 * INITIAL_STEP_SIZE = <integer>
 * The initial step size used in the ascent.
 * Default: 1
 *
 * INITIAL_TOUR_ALGORITHM = { BORUVKA | CVRP | GREEDY | MOORE | MTSP |
 *       NEAREST-NEIGHBOR | QUICK-BORUVKA | SIERPINSKI | SOP | TSPDL | WALK }
 * Specifies the algorithm for obtaining an initial tour.
 * Default: WALK
 *
 * INITIAL_TOUR_FILE = <string>
 * Specifies the name of a file containing a tour to be used as the
 * initial tour in the search. The tour is given by a list of integers
 * giving the sequence in which the nodes are visited in the tour.
 * The tour is terminated by a -1.
 * See also INITIAL_TOUR_FRACTION.
 *
 * INITIAL_TOUR_FRACTION = <real in [0;1]>
 * Specifies the fraction of the initial tour to be constructed by means
 * of INITIAL_TOUR_FILE edges.
 * Default: 1.0
 *
 * INPUT_TOUR_FILE = <string>
 * Specifies the name of a file containing a tour. The tour is used to
 * limit the search (the last edge to be excluded in a non-gainful move
 * must not belong to the tour). In addition, the Alpha field of its
 * edges is set to zero. The tour is given by a list of integers giving
 * the sequence in which the nodes are visited in the tour. The tour is
 * terminated by a -1.
 *
 * K = <integer>
 * The k-value used for solving a k-TSP instance.
 *
 * KICK_TYPE = <integer>
 * Specifies the value of k for a random k-swap kick (an extension of the
 * double-bridge move). If KICK_TYPE is zero, then the LKH's special kicking
 * strategy, WALK, is used.
 * Default: 0
 *
 * KICKS = <integer>
 * Specifies the number of times to "kick" a tour found by Lin-Kernighan.
 * Each kick is a random k-swap kick-move. However, if KICKS is zero, then
 * LKH's special kicking strategy, WALK, is used.
 * Default: 1
 *
 * MAKESPAN = { YES | NO }
 * Specifies if makespan optimization is to be used for a TSPTW instance.
 * Default: NO
 *
 * MAX_BREADTH = <integer>
 * The maximum number of candidate edges considered at each level of
 * the search for a move.
 * Default: INT_MAX
 *
 * MAX_CANDIDATES = <integer> [ SYMMETRIC ]
 * The maximum number of candidate edges to be associated with each node.
 * The integer may be followed by the keyword SYMMETRIC, signifying
 * that the candidate set is to be complemented such that every candidate
 * edge is associated with both its two end nodes.
 * If MAX_CANDIDATES is zero the candidate sets are made up of the
 * edges represented in the CANDIDATE_FILEs, the INITIAL_TOUR_FILE,
 * the INPUT_TOUR_FILE, the SUBPROBLEM_TOUR_FILE, and the MERGE_TOUR_FILEs.
 * Default: 5
 *
 * MAX_SWAPS = <integer>
 * Specifies the maximum number of swaps (flips) allowed in any search
 * for a tour improvement.
 * Default: DIMENSION
 *
 * MAX_TRIALS = <integer>
 * The maximum number of trials in each run.
 * Default: DIMENSION
 *
 * MERGE_TOUR_FILE = <string>
 * Specifies the name of a tour to be merged. The edges of the tour are
 * added to the candidate sets.
 * It is possible to give more than two MERGE_TOUR_FILE specifications.
 *
 * MOVE_TYPE = <integer> [ SPECIAL ]
 * Specifies the move type to be used in as submove in Lin-Kernighan.
 * An integer value k >= 2 signifies that a sequential k-opt move is to be used.
 * The specifier SPECIAL can be given in order to use LKH-3’s specially
 * designed moves. For this type of moves, k must be 3 or 5.
 * Default: 5
 *
 * MTSP_MIN_SIZE = <integer>
 * Specifies the minimum number of cities each salesman must visit in an
 * MTSP or CVRP instance. If negative, its value is set to
 *     DIMENSION/(ceil(1.0*DIMENSION/MTSP_MAX_SIZE) + 1)
 * Default: 1
 *
 * MTSP_MAX_SIZE = <integer>
 * Specifies the maximum number of cities each salesman may visit in an
 * MTSP or CVRP instance.
 * Default: DIMENSION-1.
 *
 * MTSP_OBJECTIVE = [ MINMAX | MINMAX_SIZE | MINSUM ]
 * Specifies the objective function type for a multiple traveling salesman
 * problem.
 *      MINMAX      - Minimize the length of the longest route.
 *      MINMAX_SIZE - Minimize the size of the largest route.
 *      MINSUM      - Minimize the length of the overall tour.
 *                    All routes must satisfy the MTSP_MIN_SIZE and
 *                    MTSP_MAX_SIZE constraints.
 *
 * MTSP_SOLUTION_FILE = <string>
 * Specifies the name of a file where the solution of an MTSP (or VRP) instance
 * is to be written.
 * The character $ in the name has a spe¬cial meaning. All occurrences are
 * replaced by the cost of the tour.
 *
 * NONSEQUENTIAL_MOVE_TYPE = <integer>
 * Specifies the nonsequential move type to be used. A value K >= 4
 * signifies that attempts are made to improve a tour by nonsequential
 * k-opt moves where 4 <= k <= K. Note, however, that the effect depends
 * on the specifications of PATCHING_C and PATCHING_A.
 * Default: (MOVE_TYPE + PATCHING_C + PATCHING_A - 1)
 *
 * OUTPUT_TOUR_FILE = <string>
 * Specifies the name of a file where the best tour is to be written.
 * Each time a trial has produced a new best tour, the tour is written
 * to this file.
 * The character $ in the name has a special meaning. All occurrences
 * are replaced by the cost of the tour.
 *
 * OPTIMUM = <integer>
 * Known optimal tour length. If STOP_AT_OPTIMUM is YES, a run will be
 * terminated if the tour length becomes equal to this value.
 * Default: MINUS_INFINITY
 *
 * PATCHING_A = <integer> [ RESTRICTED | EXTENDED ]
 * The maximum number of disjoint alternating cycles to be used for
 * patching. An attempt to patch cycles is made if the corresponding
 * non-sequential move is gainful.
 * The integer may be followed by the keyword RESTRICTED or EXTENDED.
 * The keyword RESTRICTED signifies that gainful moves are only
 * considered if all its inclusion edges are candidate edges.
 * The keyword EXTENDED signifies that the non-sequential move need
 * not be gainful if only all its inclusion edges are candidate edges.
 * Default: 1
 *
 * PATCHING_C = <integer> [ RESTRICTED | EXTENDED ]
 * The maximum number of disjoint cycles to be patched in an attempt
 * to find a feasible and gainful move. An attempt to patch cycles is
 * made if the corresponding non-sequential move is gainful.
 * The integer may be followed by the keyword RESTRICTED or EXTENDED.
 * The keyword RESTRICTED signifies that gainful moves are only
 * considered if all its inclusion edges are candidate edges.
 * The keyword EXTENDED signifies that the non-sequential move need
 * not be gainful if only all its inclusion edges are candidate edges.
 * Default: 0
 *
 * PI_FILE = <string>
 * Specifies the name of a file to which penalties (Pi-values determined
 * by the ascent) are to be written. If the file already exists, the
 * penalties are read from the file, and the ascent is skipped.
 * The first line of the file contains the number of nodes. Each of the
 * following lines is of the form
 *       <integer> <integer>
 * where the first integer is a node number, and the second integer is
 * the Pi-value associated with the node.
 * The file name "0" represents a file with all Pi-values equal to zero.
 *
 * POPMUSIC_INITIAL_TOUR = { YES | NO }
 * Specifies whether the best POPMUSIC tour is to be used as intial tour
 * for Lin-Kernighan.
 * Default: NO
 *
 * POPMUSIC_MAX_NEIGHBORS = <int>
 * Maximum number of nearest neighbors used as candidates in 3-opt for
 * POPMUSIC.
 * Default: 5
 * 
 * POPMUSIC_SAMPLE_SIZE = <int>
 * Sample size.
 * Default: 10
 *
 * POPMUSIC_SOLUTIONS = <int>
 * Number of solutions to be generated.
 * Default: 50
 *
 * POPMUSIC_TRIALS = <int>
 * Number of trials used in iterated 3-opt for POPMUSIC.
 * If the value is zero, the number of trials is the size of the subpath
 * to be optimized.
 * Default: 1
 *
 * POPULATION_SIZE = <integer>
 * Specifies the maximum size of the population in the genetic algorithm.
 * Default: 0
 *
 * PRECISION = <integer>
 * The internal precision in the representation of transformed distances:
 *    d[i][j] = PRECISION*c[i][j] + pi[i] + pi[j],
 * where d[i][j], c[i][j], pi[i] and pi[j] are all integral.
 * Default: 100 (which corresponds to 2 decimal places)
 *
 * PROBABILITY = <integer>
 * The probalilty for PTSP (in percent).
 * Default: 100.
 *
 * RECOMBINATION = { IPT | GPX2 | CLARIST }
 * Default: IPT
 *
 * RESTRICTED_SEARCH = { YES | NO }
 * Specifies whether the following search pruning technique is used:
 * The first edge to be broken in a move must not belong to the currently
 * best solution tour. When no solution tour is known, it must not belong
 * to the minimum spanning 1-tree.
 * Default: YES
 *
 * RUNS = <integer>
 * The total number of runs.
 * Default: 10
 *
 * SALESMEN = <integer>
 * Specifies the number of salesmen/vehicles.
 * Default: 1
 *
 * SCALE = <integer>
 * Scale factor. Distances are multiplied by this factor.
 * Default: 1
 *
 * SEED = <integer>
 * Specifies the initial seed for random number generation. If zero, the
 * seed is derived from the system clock.
 * Default: 1
 *
 * SINTEF_SOLUTION_FILE = <string>
 * Specifies the name of a file where the solution of an MTSP or VRP instance
 * is to be written. The solution is written in SINTEF format.
 * The character $ in the name has a special meaning. All occurrences are
 * replaced by the cost of the tour.
 *
 * SPECIAL
 * Equivalent to the following parameter specifications:
 *
 *    GAIN23 = NO
 *    KICKS = 1
 *    KICK_TYPE = 4
 *    MAX_SWAPS = 0
 *    MOVE_TYPE = 5 SPECIAL
 *    POPULATION_SIZE = 10
 *
 * STOP_AT_OPTIMUM = { YES | NO }
 * Specifies whether a run is stopped, if the tour length becomes equal
 * to OPTIMUM.
 * Default: YES
 *
 * SUBGRADIENT = { YES | NO }
 * Specifies whether the Pi-values should be determined by subgradient
 * optimization.
 * Default: YES
 *
 * SUBPROBLEM_SIZE = <integer> [ DELAUNAY | KARP | K-CENTER | K-MEANS | MOORE |
 *                               ROHE | SIERPINSKI ] [ BORDERS ] [ COMPRESSED ]
 * The number of nodes in a division of the original problem into subproblems.
 * The division is made according to the tour given by SUBPROBLEM_TOUR_FILE.
 * The value 0 signifies that no division is made.
 * By default, the subproblems are determined by subdividing the tour into
 * segments of equal size. However, the integer may be followed by DELAUNAY,
 * KARP, K-CENTER, K-MEANS, MOORE, ROHE or SIERPINSKI. DELAUNAY specifies that
 * the Delaunay partitioning scheme is used, KARP that Karp's partitioning
 * scheme is used, K-CENTER that a partitioning scheme based on K-center
 * clustering, K-MEANS that a partitioning scheme based on K-means clustering
 * is used, ROHE that Rohe's random rectangle/cube partitioning scheme is used,
 * and MOORE or SIERPINSKI that a partitioning scheme based on either a Moore
 * or Sierpinski space-filling curve is used.
 * The BORDERS specification signifies that the subproblems along the borders
 * between subproblems are to be solved too.
 * The COMPRESSED specification signifies that each subproblem is compressed by
 * removing from the problem all nodes with two incident subproblem tour edges
 * that belong to all tours to be merged (at least two MERGE_TOUR_FILEs should
 * be given).
 * Default: 0
 *
 * SUBPROBLEM_TOUR_FILE = <string>
 * Specifies the name of a file containing a tour to be used for dividing
 * the original problem into subproblems. The approximate number of nodes
 * in each is * given by SUBPROBLEM_SIZE.
 * The tour is given by a list of integers giving the sequence in which the
 * nodes are visited in the tour. The tour is terminated by a -1
 *
 * SUBSEQUENT_MOVE_TYPE = <integer> [ SPECIAL ]
 * Specifies the move type to be used for all moves following the first move
 * in a sequence of moves. The value K >= 2 signifies that a K-opt move is to
 * be used. The value 0 signifies that all moves are of the same type
 * (K = MOVE_TYPE).
 * The specifier SPECIAL can be given in order to use LKH-3’s specially
 * designed moves. For this type of moves, k must be 3 or 5.
 * Default: 0
 *
 * SUBSEQUENT_PATCHING = { YES | NO }
 * Specifies whether patching is used for moves following the first move
 * in a sequence of moves.
 * Default: YES
 *
 * TIME_LIMIT = <real>
 * Specifies a time limit in seconds for each run.
 * Default: DBL_MAX
 *
 * TOTAL_TIME_LIMIT = <real>
 * Specifies a total time limit in seconds.
 * Default: DBL_MAX
 *
 * TOUR_FILE = <string>
 * Specifies the name of a file where the best tour is to be written.
 * When a run has produced a new best tour, the tour is written to this file.
 * The character $ in the name has a special meaning. All occurrences
 * are replaced by the cost of the tour.
 *
 * TRACE_LEVEL = <integer>
 * Specifies the level of detail of the output given during the solution
 * process. The value 0 signifies a minimum amount of output. The higher
 * the value is the more information is given.
 * Default: 1
 *
 * VEHICLES = <integer>
 * Specifies the number of vehicles/salesmen.
 * Default: 1
 * 
 * List of abbreviations
 * ---------------------
 *
 * A string value may be abbreviated to the first few letters of the string,
 * if that abbreviation is unambiguous.
 *
 *     Value        Abbreviation
 *     ALPHA             A
 *     BORDERS           B
 *     BORUVKA           B
 *     CLARIST           C
 *     COMPRESSED        C
 *     CVRP              C
 *     DELAUNAY          D
 *     EXTENDED          E
 *     GPX2              G
 *     GREEDY            G
 *     IPT               I
 *     KARP              KA
 *     K-CENTER          K-C
 *     K-MEANS           K-M
 *     MOORE             MO
 *     MTSP              MT
 *     NEAREST-NEIGHBOR  N
 *     NO                N
 *     POPMUSIC          P
 *     PURE              P
 *     QUADRANT          Q
 *     QUICK-BORUVKA     Q
 *     RESTRICTED        R
 *     ROHE              R
 *     SIERPINSKI        SI
 *     SOP               SO
 *     SPECIAL           S
 *     SYMMETRIC         S
 *     TSPDL             T
 *     WALK              W
 *     YES               Y
 */

    /* 可设置的参数 */
    ProblemFileName = PiFileName = InputTourFileName =
    OutputTourFileName = TourFileName = 0;
    CandidateFiles = MergeTourFiles = 0;
    AscentCandidates = 50;
    BackboneTrials = 0;
    Backtracking = 0;
    BWTSP_B = 0;
    BWTSP_Q = 0;
    BWTSP_L = INT_MAX;
    CandidateSetSymmetric = 0;
    // CandidateSetType = ALPHA;
    CandidateSetType = POPMUSIC;
    Crossover = ERXT;
    DelaunayPartitioning = 0;
    DelaunayPure = 0;
    DemandDimension = 1;
    DistanceLimit = DBL_MAX;
    Excess = -1;
    ExternalSalesmen = 0;
    ExtraCandidates = 0;
    ExtraCandidateSetSymmetric = 0;
    ExtraCandidateSetType = QUADRANT;
    Gain23Used = 1;
    GainCriterionUsed = 1;
    GridSize = 1000000.0;
    InitialPeriod = -1;
    InitialStepSize = 0;
    InitialTourAlgorithm = WALK;
    InitialTourFraction = 1.0;
    KarpPartitioning = 0;
    KCenterPartitioning = 0;
    KMeansPartitioning = 0;
    Kicks = 1;
    KickType = 0;
    MaxBreadth = INT_MAX;
    MaxCandidates = 5;
    MaxPopulationSize = 0;
    MaxSwaps = -1;
    MaxTrials = -1;
    MoorePartitioning = 0;
    MoveType = 5;
    MoveTypeSpecial = 0;
    MTSPDepot = 1;
    MTSPMinSize = 1;
    MTSPMaxSize = -1;
    MTSPObjective = -1;
    NonsequentialMoveType = -1;
    Optimum = MINUS_INFINITY;
    PatchingA = 1;
    PatchingC = 0;
    PatchingAExtended = 0;
    PatchingARestricted = 0;
    PatchingCExtended = 0;
    PatchingCRestricted = 0;
    Precision = 100;
    Probability = 100;
    // POPMUSIC_InitialTour = 0;
    POPMUSIC_InitialTour = 1;
    POPMUSIC_MaxNeighbors = 5;
    POPMUSIC_SampleSize = 10;
    POPMUSIC_Solutions = 50;
    // POPMUSIC_Solutions = 25; // 这个减半后，IO为10000的也能20秒ascent完，且因为没超时，最终分数还会高些
    POPMUSIC_Trials = 1;
    Recombination = IPT;
    RestrictedSearch = 1;
    RohePartitioning = 0;
    Runs = 0;
    Salesmen = 1;
    Scale = -1;
    Seed = 1;
    SierpinskiPartitioning = 0;
    StopAtOptimum = 1;
    // Subgradient = 1;
    Subgradient = 0;
    SubproblemBorders = 0;
    SubproblemsCompressed = 0;
    SubproblemSize = 0;
    SubsequentMoveType = 0;
    SubsequentMoveTypeSpecial = 0;
    SubsequentPatching = 1;
    TimeLimit = DBL_MAX;
    TotalTimeLimit = 20;
    TraceLevel = 1;
    TSPTW_Makespan = 0;

    /* other */
    MaxMatrixDimension = 20006;
    MergeWithTour =
        Recombination == GPX2 ? MergeWithTourGPX2 :
        Recombination == CLARIST ? MergeWithTourCLARIST :
                                   MergeWithTourIPT;
}

void ReadProblem(int **adjMat, int MatDimension)
{   
   FreeStructures();
    FirstNode = 0;
    C = 0;
    c = 0;

    /* (1) The specification part 定义问题类型、规模和数据输入格式 */
    Name = Copy("IOSchedule");
    Type = Copy("ATSP");
    ProblemType = ATSP;
    Asymmetric = 1;
    WeightType = EXPLICIT;
    WeightFormat = FULL_MATRIX;
    EdgeWeightType = Copy("EXPLICIT");
    EdgeWeightFormat = Copy("FULL_MATRIX");
    EdgeDataFormat = NULL;
    NodeCoordType = Copy("NO_COORDS");
    DisplayDataType = 0;
    Distance = Distance_EXPLICIT;;

    Dimension = MatDimension;
    DimensionSaved = Dim = Dimension;
    DistanceLimit = DBL_MAX;
    Salesmen = 1;

    /* (2) The data part 根据上方定义，进行实际数据输入 */
    /* NODE_COORD_SECTION */
    // NULL
    /* EDGE_DATA_SECTION */
    // NULL
    /* FIXED_EDGES_SECTION */
    // NULL
    /* DISPLAY_DATA_SECTION */
    // NULL
    /* TOUR_SECTION */
    // NULL
    /* EDGE_WEIGHT_SECTION */
    Read_EDGE_WEIGHT_SECTION(adjMat);

    Swaps = 0;

    /* Adjust parameters 调整参数。比如某个参数的默认值依赖于问题类型或输入数据，
       在此前(ReadParameter 函数中)无法确定值是多少，仅是设定为一个默认值标记，
       所以在此处(ReadProblem 函数最后)其默认值才被实际的设定 */
    int i, j, K;
    if (Seed == 0)
        Seed = (unsigned)time(0);
    if (InitialStepSize == 0)
        InitialStepSize = 1;
    if (MaxSwaps < 0)
        MaxSwaps = Dimension;
    if (KickType > Dimension / 2)
        KickType = Dimension / 2;
    if (Runs == 0)
        Runs = 10;
    if (MaxCandidates > Dimension - 1)
        MaxCandidates = Dimension - 1;
    if (ExtraCandidates > Dimension - 1)
        ExtraCandidates = Dimension - 1;
    if (SubproblemSize >= Dimension)
        SubproblemSize = Dimension;
    else if (SubproblemSize == 0) {
        if (AscentCandidates > Dimension - 1)
            AscentCandidates = Dimension - 1;
        if (InitialPeriod < 0) {
            InitialPeriod = Dimension / 2;
            if (InitialPeriod < 100)
                InitialPeriod = 100;
        }
        if (Excess < 0)
            Excess = 1.0 / DimensionSaved * Salesmen;
        if (MaxTrials == -1)
            MaxTrials = Dimension;
        HeapMake(Dimension);
    }
    if (POPMUSIC_MaxNeighbors > Dimension - 1)
        POPMUSIC_MaxNeighbors = Dimension - 1;
    if (POPMUSIC_SampleSize > Dimension)
        POPMUSIC_SampleSize = Dimension;
    Depot = &NodeSet[MTSPDepot];
    if (ProblemType == CVRP) {
        Node *N;
        int MinSalesmen;
        if (Capacity <= 0)
            eprintf("CAPACITY not specified");
        TotalDemand = 0;
        N = FirstNode;
        do
            TotalDemand += N->Demand;
        while ((N = N->Suc) != FirstNode);
        MinSalesmen =
            TotalDemand / Capacity + (TotalDemand % Capacity != 0);
        if (Salesmen == 1) {
            Salesmen = MinSalesmen;
            if (Salesmen > Dimension)
                eprintf("CVRP: SALESMEN larger than DIMENSION");
        } else if (Salesmen < MinSalesmen)
            eprintf("CVRP: SALESMEN too small to meet demand");
        assert(Salesmen >= 1 && Salesmen <= Dimension);
        if (Salesmen == 1)
            ProblemType = TSP;
        Penalty = Penalty_CVRP;
    } else if (ProblemType == SOP || ProblemType == M1_PDTSP) {
        Constraint *Con;
        Node *Ni, *Nj;
        int n, k;
        OldDistance = Distance;
        Distance = Distance_SOP;
        if (ProblemType == M1_PDTSP) {
            for (i = 2; i < Dim; i++) {
                Ni = &NodeSet[i];
                for (k = n = 0; k < DemandDimension; k++) {
                    n = Ni->M_Demand[k];
                    if (n >= 0)
                        continue;
                    for (j = 2; j < Dim; j++) {
                        if (j == i)
                            continue;
                        Nj = &NodeSet[j];
                        if (Nj->M_Demand[k] == -n) {
                            Ni->C[j] = -1;
                            break;
                        }
                    }
                }
            }
        }
        for (j = 2; j < Dim; j++) {
            Nj = &NodeSet[j];
            for (i = 2; i < Dim; i++) {
                if (i != j && Nj->C[i] == -1) {
                    Ni = &NodeSet[i];
                    assert(Con =
                        (Constraint *)malloc(sizeof(Constraint)));
                    Con->t1 = Ni;
                    Con->t2 = Nj;
                    Con->Suc = FirstConstraint;
                    FirstConstraint = Con;
                    Con->Next = Ni->FirstConstraint;
                    Ni->FirstConstraint = Con;
                }
            }
        }
        Salesmen = 1;
        Penalty = ProblemType == SOP ? Penalty_SOP : Penalty_M1_PDTSP;
    }
    if (ProblemType == TSPTW) {
        Salesmen = 1;
        Penalty = Penalty_TSPTW;
    } else
        TSPTW_Makespan = 0;
    if (Salesmen > 1) {
        if (Salesmen > Dim)
            eprintf("Too many salesmen/vehicles (> DIMENSION)");
        MTSP2TSP();
    }
    if (ExternalSalesmen > Salesmen)
        ExternalSalesmen = Salesmen;
    if (ProblemType == ACVRP)
        Penalty = Penalty_ACVRP;
    else if (ProblemType == CCVRP)
        Penalty = Penalty_CCVRP;
    else if (ProblemType == CTSP)
        Penalty = Penalty_CTSP;
    else if (ProblemType == CVRPTW)
        Penalty = Penalty_CVRPTW;
    else if (ProblemType == MLP)
        Penalty = Penalty_MLP;
    else if (ProblemType == OVRP)
        Penalty = Penalty_OVRP;
    else if (ProblemType == PDTSP)
        Penalty = Penalty_PDTSP;
    else if (ProblemType == PDTSPF)
        Penalty = Penalty_PDTSPF;
    else if (ProblemType == PDTSPL)
        Penalty = Penalty_PDTSPL;
    else if (ProblemType == PDPTW)
        Penalty = Penalty_PDPTW;
    else if (ProblemType == ONE_PDTSP)
        Penalty = Penalty_1_PDTSP;
    else if (ProblemType == M_PDTSP)
        Penalty = Penalty_M_PDTSP;
    else if (ProblemType == M1_PDTSP)
        Penalty = Penalty_M1_PDTSP;
    else if (ProblemType == RCTVRP || ProblemType == RCTVRPTW)
        Penalty = Penalty_RCTVRP;
    else if (ProblemType == TRP)
        Penalty = Penalty_TRP;
    else if (ProblemType == TSPDL)
        Penalty = Penalty_TSPDL;
    else if (ProblemType == TSPPD)
        Penalty = Penalty_TSPPD;
    if (ProblemType == VRPB)
        Penalty = Penalty_VRPB;
    else if (ProblemType == VRPBTW)
        Penalty = Penalty_VRPBTW;
    else if (ProblemType == VRPPD)
        Penalty = Penalty_VRPPD;
    if (BWTSP_B > 0) {
        if (Penalty)
            eprintf("BWTSP not compatible with problem type %s\n", Type);
        ProblemType = BWTSP;
        free(Type);
        Type = Copy("BWTSP");
        Penalty = Penalty_BWTSP;
        if (BWTSP_L != INT_MAX)
            BWTSP_L *= Scale;
    }
    if (Penalty && (SubproblemSize > 0 || SubproblemTourFile))
        eprintf("Partitioning not implemented for constrained problems");
    Depot->DepotId = 1;
    for (i = Dim + 1; i <= DimensionSaved; i++)
        NodeSet[i].DepotId = i - Dim + 1;
    if (Dimension != DimensionSaved) {
        NodeSet[Depot->Id + DimensionSaved].DepotId = 1;
        for (i = Dim + 1; i <= DimensionSaved; i++)
            NodeSet[i + DimensionSaved].DepotId = i - Dim + 1;
    }
    if (ServiceTime != 0) {
        for (i = 1; i <= Dim; i++)
            NodeSet[i].ServiceTime = ServiceTime;
        Depot->ServiceTime = 0;
    }
    if (CostMatrix == 0 && Dimension <= MaxMatrixDimension &&
        Distance != 0 && Distance != Distance_1 && Distance != Distance_LARGE
        && Distance != Distance_LARGE && Distance != Distance_ATSP
        && Distance != Distance_MTSP && Distance != Distance_SPECIAL) {
        Node *Ni, *Nj;
        assert(CostMatrix =
            (int *)calloc((size_t)Dim * (Dim - 1) / 2, sizeof(int)));
        Ni = FirstNode->Suc;
        do {
            Ni->C =
                &CostMatrix[(size_t)(Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
            if (ProblemType != HPP || Ni->Id <= Dim)
                for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                    Ni->C[Nj->Id] = Fixed(Ni, Nj) ? 0 : Distance(Ni, Nj);
            else
                for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                    Ni->C[Nj->Id] = 0;
        } while ((Ni = Ni->Suc) != FirstNode);
        c = 0;
        WeightType = EXPLICIT;
    }
    if (ProblemType == TSPTW ||
        ProblemType == CVRPTW || ProblemType == VRPBTW ||
        ProblemType == PDPTW || ProblemType == RCTVRPTW) {
        M = INT_MAX / 2 / Precision;
        for (i = 1; i <= Dim; i++) {
            Node *Ni = &NodeSet[i];
            for (j = 1; j <= Dim; j++) {
                Node *Nj = &NodeSet[j];
                if (Ni != Nj &&
                    Ni->Earliest + Ni->ServiceTime + Ni->C[j] > Nj->Latest)
                    Ni->C[j] = M;
            }
        }
    }
    C = WeightType == EXPLICIT ? C_EXPLICIT : C_FUNCTION;
    D = WeightType == EXPLICIT ? D_EXPLICIT : D_FUNCTION;
    if (ProblemType != CVRP && ProblemType != CVRPTW &&
        ProblemType != CTSP &&
        ProblemType != TSP && ProblemType != ATSP) {
        M = INT_MAX / 2 / Precision;
        for (i = Dim + 1; i <= DimensionSaved; i++) {
            for (j = 1; j <= DimensionSaved; j++) {
                if (j == i)
                    continue;
                if (j == MTSPDepot || j > Dim)
                    NodeSet[i].C[j] = NodeSet[MTSPDepot].C[j] = M;
                NodeSet[i].C[j] = NodeSet[MTSPDepot].C[j];
                NodeSet[j].C[i] = NodeSet[j].C[MTSPDepot];
            }
        }
        if (ProblemType == CCVRP || ProblemType == OVRP)
            for (i = 1; i <= Dim; i++)
                NodeSet[i].C[MTSPDepot] = 0;
    }
    if (Precision > 1 && CostMatrix) {
        for (i = 2; i <= Dim; i++) {
            Node *N = &NodeSet[i];
            for (j = 1; j < i; j++)
                if (N->C[j] * Precision / Precision != N->C[j])
                    eprintf("PRECISION (= %d) is too large", Precision);
        }
    }
    if (SubsequentMoveType == 0) {
        SubsequentMoveType = MoveType;
        SubsequentMoveTypeSpecial = MoveTypeSpecial;
    }
    K = MoveType >= SubsequentMoveType || !SubsequentPatching ?
        MoveType : SubsequentMoveType;
    if (PatchingC > K)
        PatchingC = K;
    if (PatchingA > 1 && PatchingA >= PatchingC)
        PatchingA = PatchingC > 2 ? PatchingC - 1 : 1;
    if (NonsequentialMoveType == -1 ||
        NonsequentialMoveType > K + PatchingC + PatchingA - 1)
        NonsequentialMoveType = K + PatchingC + PatchingA - 1;
    if (PatchingC >= 1 && NonsequentialMoveType >= 4) {
        BestMove = BestSubsequentMove = BestKOptMove;
        if (!SubsequentPatching && SubsequentMoveType <= 5) {
            MoveFunction BestOptMove[] =
            { 0, 0, Best2OptMove, Best3OptMove,
                Best4OptMove, Best5OptMove
            };
            BestSubsequentMove = BestOptMove[SubsequentMoveType];
        }
    } else {
        MoveFunction BestOptMove[] = { 0, 0, Best2OptMove, Best3OptMove,
            Best4OptMove, Best5OptMove
        };
        BestMove = MoveType <= 5 ? BestOptMove[MoveType] : BestKOptMove;
        BestSubsequentMove = SubsequentMoveType <= 5 ?
            BestOptMove[SubsequentMoveType] : BestKOptMove;
    }
    if (MoveTypeSpecial)
        BestMove = BestSpecialOptMove;
    if (SubsequentMoveTypeSpecial)
        BestSubsequentMove = BestSpecialOptMove;
    if (ProblemType == HCP || ProblemType == HPP)
        MaxCandidates = 0;
    if (TraceLevel >= 1) {
        printff("done\n");
        PrintParameters();
    }
    free(LastLine);
    LastLine = 0;
}

void OutputTourResult(int* tourResult, int *tourCost){
    int i, j, n, Forwards;

    n = ProblemType != ATSP ? Dimension : Dimension / 2;
    for (i = 1; i < n && BestTour[i] != 1; i++);
    Forwards = ProblemType == ATSP ||
        BestTour[i < n ? i + 1 : 1] < BestTour[i > 1 ? i - 1 : Dimension];
    for (j = 1; j <= n; j++) {
        tourResult[j - 1] = BestTour[i];
        if (Forwards) {
            if (++i > n)
                i = 1;
        } else if (--i < 1)
            i = n;
    }
    *tourCost = BestCost;
}

void ReSetLKH(){
    // [TODO 1: 释放 LKH内核 申请(malloc等)的内存, 可以对照源文件 AllocateStructures.c 来做]
    /* tips: 如果在程序结束前没有调用 free，申请的内存将在程序退出时被操作系统回收，
        但这并不算是良好的编程习惯。确保释放内存是管理动态内存的好方法。 */

    // [TODO 2: 重置 LKH内核 的全局变量, 可以对照源文件 LKH.c 来做]
    /* tips: 因为 LKH内核 的运行依赖于全局变量，如果在一个进程内存在多次调用 LKH内核 的话，
        除第一次调用外，后续调用需要重置全局变量，这样才能确保不会对后续调用产生影响 */

    // [TODO 3: 重置 LKH内核 涉及的静态变量, 这些静态变量分布于多个源文件中，需要逐个确认]
    /* tips: 因为 LKH内核 的运行同样依赖于这些静态变量，如果在一个进程内存在多次调用 LKH内核 的话，
        除第一次调用外，后续调用需要重置这些静态变量，这样才能确保不会对后续调用产生影响 */
    /* method: 由于静态变量的作用域仅限于其定义的源文件，无法在其他源文件中直接重新初始化。
        因此，需要使用函数接口：即在定义静态变量的源文件内部，创建一个函数来重新初始化其涉及的静态变量，
        其他源文件就可以通过调用这个函数来实现重新初始化。 */
}

int solveTSP(int **adjMat, int MatDimension, int *tourResult, int *tourCost){
    /* 参数设定 */
    ReadParameters();
    /* 问题信息设定、数据读入 */
    ReadProblem(adjMat, MatDimension);
    /* 调用 LKH内核 */
    int ret = LKHmain();
    /* 提取求解结果 */
    OutputTourResult(tourResult, tourCost);
    /* 重置 LKH内核 状态*/
    ReSetLKH();
    return ret;
}

static char *Copy(char *S)
{
    char *Buffer;

    if (!S || strlen(S) == 0)
        return 0;
    Buffer = (char *) malloc(strlen(S) + 1);
    strcpy(Buffer, S);
    return Buffer;
}
static void CheckSpecificationPart()
{
    if (ProblemType == -1)
        eprintf("TYPE is missing");
    if (Dimension < 3)
        eprintf("DIMENSION < 3 or not specified");
    if (WeightType == -1 && !Asymmetric && ProblemType != HCP &&
        ProblemType != HPP && !EdgeWeightType && ProblemType != STTSP)
        eprintf("EDGE_WEIGHT_TYPE is missing");
    if (WeightType == EXPLICIT && WeightFormat == -1 && !EdgeWeightFormat)
        eprintf("EDGE_WEIGHT_FORMAT is missing");
    if (WeightType == EXPLICIT && WeightFormat == FUNCTION)
        eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if (WeightType != EXPLICIT &&
        (WeightType != SPECIAL || CoordType != NO_COORDS) &&
        WeightType != -1 && WeightFormat != -1 && WeightFormat != FUNCTION)
        eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if ((ProblemType == ATSP || ProblemType == SOP) &&
        WeightType != EXPLICIT && WeightType != -1)
        eprintf("Conflicting TYPE and EDGE_WEIGHT_TYPE");
    if (CandidateSetType == DELAUNAY && !TwoDWeightType() &&
        MaxCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = DELAUNAY");
    if (CandidateSetType == NN && !TwoDWeightType() &&
        MaxCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for "
             "CANDIDATE_SET_TYPE = NEAREST-NEIGHBOR");
    if (CandidateSetType == QUADRANT && !TwoDWeightType() &&
        !ThreeDWeightType() && MaxCandidates + ExtraCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = QUADRANT");
    if (ExtraCandidateSetType == QUADRANT && !TwoDWeightType() &&
        !ThreeDWeightType() && ExtraCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for EXTRA_CANDIDATE_SET_TYPE = "
             "QUADRANT");
    if (InitialTourAlgorithm == QUICK_BORUVKA && !TwoDWeightType() &&
        !ThreeDWeightType())
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for INITIAL_TOUR_ALGORITHM = "
             "QUICK-BORUVKA");
    if (InitialTourAlgorithm == SIERPINSKI && !TwoDWeightType())
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for INITIAL_TOUR_ALGORITHM = "
             "SIERPINSKI");
    if (DelaunayPartitioning && !TwoDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for DELAUNAY specification");
    if (KarpPartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for KARP specification");
    if (KCenterPartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for K-CENTER specification");
    if (KMeansPartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for K-MEANS specification");
    if (MoorePartitioning && !TwoDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for MOORE specification");
    if (RohePartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for ROHE specification");
    if (SierpinskiPartitioning && !TwoDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for SIERPINSKI specification");
    if (SubproblemBorders && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for BORDERS specification");
    if (InitialTourAlgorithm == MTSP_ALG && Asymmetric)
        eprintf("INTIAL_TOUR_ALGORITHM = MTSP is not applicable for "
                "asymetric problems");
}
static void CreateNodes()
{
    Node *Prev = 0, *N = 0;
    int i;
    if (Dimension <= 0)
        eprintf("DIMENSION is not positive (or not specified)");
    if (Asymmetric) {
        Dim = DimensionSaved;
        DimensionSaved = Dimension + Salesmen - 1;
        Dimension = 2 * DimensionSaved;
    } else if (ProblemType == HPP) {
        Dimension++;
        if (Dimension > MaxMatrixDimension)
            eprintf("DIMENSION too large in HPP problem");
    }
    NodeSet = (Node *) calloc(Dimension + 1, sizeof(Node));
    for (i = 1; i <= Dimension; i++, Prev = N) {
        N = &NodeSet[i];
        if (i == 1)
            FirstNode = N;
        else
            Link(Prev, N);
        N->Id = N->OriginalId = i;
        if (MergeTourFiles >= 1)
            N->MergeSuc = (Node **) calloc(MergeTourFiles, sizeof(Node *));
        N->Earliest = 0;
        N->Latest = INT_MAX;
    }
    Link(N, FirstNode);
}
static int FixEdge(Node * Na, Node * Nb)
{
    if (!Na->FixedTo1 || Na->FixedTo1 == Nb)
        Na->FixedTo1 = Nb;
    else if (!Na->FixedTo2 || Na->FixedTo2 == Nb)
        Na->FixedTo2 = Nb;
    else
        return 0;
    if (!Nb->FixedTo1 || Nb->FixedTo1 == Na)
        Nb->FixedTo1 = Na;
    else if (!Nb->FixedTo2 || Nb->FixedTo1 == Na)
        Nb->FixedTo2 = Na;
    else
        return 0;
    return 1;
}
static void Read_EDGE_WEIGHT_SECTION(int **adjMat)
{
    Node *Ni;
    int i, j, n, W;
    double w;
    CheckSpecificationPart();
    if (!FirstNode)
        CreateNodes(); // 创建结点数组，是一个头尾相连的双链表数据结构，FirstNode指向第一个结点
    if (!Asymmetric) {
        CostMatrix = (int *) calloc((size_t) Dimension * (Dimension - 1) / 2,
                                    sizeof(int));
        Ni = FirstNode->Suc;
        do {
            Ni->C =
                &CostMatrix[(size_t) (Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
        }
        while ((Ni = Ni->Suc) != FirstNode);
    } else {
        n = Dimension / 2;
        CostMatrix = (int *) calloc((size_t) n * n, sizeof(int));
        for (Ni = FirstNode; Ni->Id <= n; Ni = Ni->Suc)
            Ni->C = &CostMatrix[(size_t) (Ni->Id - 1) * n] - 1;
    }
    if (ProblemType == HPP)
        Dimension--;
    if (Scale < 1)
        Scale = 1;
    // 写入矩阵中
    for (i = 1; i <= Dim; i++) {
        Ni = &NodeSet[i];
        for (j = 1; j <= Dim; j++) {
            W = adjMat[i - 1][j - 1]; // 对应位置
            if (Asymmetric) {
                Ni->C[j] = W;
                if (j != i && W > M)
                    M = W;
            } else if (j < i)
                Ni->C[j] = W;
        }
    }
    if (ProblemType == HPP)
        Dimension++;
    if (Asymmetric) {
        for (i = 1; i <= DimensionSaved; i++)
            FixEdge(&NodeSet[i], &NodeSet[i + DimensionSaved]);
        Distance = Distance_ATSP;
        WeightType = -1;
    }
}
static int TwoDWeightType()
{
    if (Asymmetric)
        return 0;
    return WeightType == EUC_2D || WeightType == MAX_2D ||
        WeightType == MAN_2D || WeightType == CEIL_2D ||
        WeightType == FLOOR_2D ||
        WeightType == GEO || WeightType == GEOM ||
        WeightType == GEO_MEEUS || WeightType == GEOM_MEEUS ||
        WeightType == ATT || WeightType == TOR_2D ||
        (WeightType == SPECIAL && CoordType == TWOD_COORDS);
}
static int ThreeDWeightType()
{
    if (Asymmetric)
        return 0;
    return WeightType == EUC_3D || WeightType == MAX_3D ||
        WeightType == MAN_3D || WeightType == CEIL_3D ||
        WeightType == FLOOR_3D || WeightType == TOR_3D ||
        WeightType == XRAY1 || WeightType == XRAY2 ||
        (WeightType == SPECIAL && CoordType == THREED_COORDS);
}
