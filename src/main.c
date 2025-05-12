#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#define SIZE 9
#define BOX 3

void fillGrid(int grid[SIZE][SIZE], bool fixed[SIZE][SIZE]);
bool fixedValueInBox(int grid[SIZE][SIZE], bool fixed[SIZE][SIZE], int blockRow, int blockCol, int val);
void initFixed(int grid[SIZE][SIZE], bool fixed[SIZE][SIZE]);
void printGrid(int grid[SIZE][SIZE]);
void generate_neighbor_inplace(int grid[SIZE][SIZE], bool fixed[SIZE][SIZE]);
int computeCost(int grid[SIZE][SIZE]);
void solveSudoku_SA(int grid[SIZE][SIZE], bool fixed[SIZE][SIZE], double T_start, double T_end, double alpha, int max_iterations);


int main(){
    int grid[SIZE][SIZE] = {
        {5, 3, 0, 0, 7, 0, 0, 0, 0},
        {6, 0, 0, 1, 9, 5, 0, 0, 0},
        {0, 9, 8, 0, 0, 0, 0, 6, 0},
        {8, 0, 0, 0, 6, 0, 0, 0, 3},
        {4, 0, 0, 8, 0, 3, 0, 0, 1},
        {7, 0, 0, 0, 2, 0, 0, 0, 6},
        {0, 6, 0, 0, 0, 0, 2, 8, 0},
        {0, 0, 0, 4, 1, 9, 0, 0, 5},
        {0, 0, 0, 0, 8, 0, 0, 7, 9}
    };


    bool fixed[SIZE][SIZE];
    initFixed(grid, fixed);

    fillGrid(grid, fixed);

    printf("initial filled grid:\n"); //informational
    printGrid(grid);
    

    // SA parameters
    double T_start = 1000.0; // temperature at the start
    double T_end = 0.01;
    double alpha = 0.9995; // 
    int max_iterations = 500000;

    solveSudoku_SA(grid, fixed, T_start, T_end, alpha, max_iterations);

    printf("\nsolved grid:\n");
    printGrid(grid);

    getchar();
    getchar();

    return 0;
}

// fills the spots; '0's on the board with numbers 1-9
// makes sure no digit is duplicated within the box
void fillGrid(int grid[SIZE][SIZE], bool fixed[SIZE][SIZE]) {

    for (int blockRow = 0; blockRow < BOX; blockRow++) {
        for (int blockCol = 0; blockCol < BOX; blockCol++) {
            int digits[SIZE];
            int count = 0;

            // collect digits 1–9
            for (int i = 0; i < SIZE; i++) {
                digits[i] = i + 1;
            }

            // shuffle digits
            for (int i = SIZE - 1; i > 0; i--) {
                int j = rand() % (i + 1);
                int temp = digits[i];
                digits[i] = digits[j];
                digits[j] = temp;
            }

            // fill unfixed cells with shuffled digits, skip over fixed
            int index = 0;
            for (int i = 0; i < BOX; i++) {
                for (int j = 0; j < BOX; j++) {
                    int row = blockRow * BOX + i;
                    int col = blockCol * BOX + j;

                    if (!fixed[row][col]) {
                        // skip used values by fixed cells in this box
                        while (index < SIZE && fixedValueInBox(grid, fixed, blockRow, blockCol, digits[index])) {
                            index++;
                        }
                        //fills the cell with the next available digit.
                        if (index < SIZE) {
                            grid[row][col] = digits[index++];
                        }
                    }
                }
            }

        }
    }
}

// check if digit already used in this box (for fixed cells)
bool fixedValueInBox(int grid[SIZE][SIZE], bool fixed[SIZE][SIZE], int blockRow, int blockCol, int val) {
    for (int i = 0; i < BOX; i++) {
        for (int j = 0; j < BOX; j++) {
            int row = blockRow * BOX + i;
            int col = blockCol * BOX + j;
            if (fixed[row][col] && grid[row][col] == val) {
                return true;
            }
        }
    }
    return false;
}



//function which initializes the bool fixed[SIZE][SIDE] board
void initFixed(int grid[SIZE][SIZE], bool fixed[SIZE][SIZE]) {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (grid[i][j] != 0)
                fixed[i][j] = true;   // 'frozen'
            else
                fixed[i][j] = false;  // can be changed
        }
    }
}


// printing the grid
void printGrid(int grid[SIZE][SIZE]) {
  for (int i = 0; i < SIZE; i++) {
  
    if (i % BOX == 0) {
    printf("--------------------------------\n");
    }
  
    for (int j = 0; j < SIZE; j++) {
  
    if (j % BOX == 0) {
        printf("| ");
    }
  
    if (grid[i][j] == 0) {
        printf(" 0 ");
    } 
    else if (grid[i][j] <= 9) {
        printf("%2d ", grid[i][j]);
    } 
    else {
        // it will replace 10-16 with a char of A-F
        printf(" %c ", 'A' + (grid[i][j] - 10));
        }
    }
      printf("\n");
    }
}


//generating neighbor

typedef struct {
    int row;
    int col;
} Cell;

// swapping two cells
void generate_neighbor_inplace(int grid[SIZE][SIZE], bool fixed[SIZE][SIZE]) {
    // chooses a random box in the grid; 
    // there are 3x3 boxes of (3x3) boxes
    int block_row = rand() % BOX;
    int block_col = rand() % BOX;

    //calculates the starting point of the box
    int start_row = block_row * BOX;
    int start_col = block_col * BOX;

    // [9][2] - 9 because the max number of possible "pairs" to switch are 9 in a box
    //          it holds a pair like (0, 1), coordinates which ones to switch
    Cell variable_cells[9];
    int count = 0;

    // collecting not frozen cells
    for (int i = 0; i < BOX; i++) {
        for (int j = 0; j < BOX; j++) {
            int row = start_row + i;
            int col = start_col + j;

            if (!fixed[row][col]) {
                variable_cells[count].row = row;
                variable_cells[count].col = col;
                count++;
            }
        }
    }

    // if fewer than 2 modifiable cells;
    if (count < 2) return;

    // pick 2 random modificable cells to swap
    int idex1 = rand() % count;
    int idex2;
    do {
        idex2 = rand() % count;
    } while (idex2 == idex1);

    // swap the values
    Cell c1 = variable_cells[idex1];
    Cell c2 = variable_cells[idex2];

    int temp = grid[c1.row][c1.col];
    grid[c1.row][c1.col] = grid[c2.row][c2.col];
    grid[c2.row][c2.col] = temp;
}


// calculating the cost = sum of "errors" a board has
int computeCost(int grid[SIZE][SIZE]) {
    int cost = 0;

    // rows
    for (int i = 0; i < SIZE; i++) {
        int count[SIZE + 1] = {0};
        for (int j = 0; j < SIZE; j++) {
            count[grid[i][j]]++;
        }
        for (int k = 1; k <= SIZE; k++) {
            if (count[k] > 1){
                cost += count[k] - 1;
            }
        }
    }

    // columns
    for (int j = 0; j < SIZE; j++) {
        int count[SIZE + 1] = {0};
        for (int i = 0; i < SIZE; i++) {
            count[grid[i][j]]++;
        }
        for (int k = 1; k <= SIZE; k++) {
            if (count[k] > 1){
                cost += count[k] - 1;
            }
        }
    }

    return cost;
}


// sudoku solver - simulated annealing
void solveSudoku_SA(int grid[SIZE][SIZE], bool fixed[SIZE][SIZE], double T_start, double T_end, double alpha, int max_iterations) {
    srand(time(NULL));

    int currentCost = computeCost(grid);
    double T = T_start;

    // LOCAL BESTS; these ones will update
    int bestCost = currentCost;
    int bestGrid[SIZE][SIZE];
    int backup[SIZE][SIZE];

    // current grid as best initial for now
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            bestGrid[i][j] = grid[i][j];
        }
    }

    int iteration = 0;
    while (T > T_end && iteration < max_iterations && bestCost > 0) {
        // backing up the current grid
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                backup[i][j] = grid[i][j];
            }
        }

        generate_neighbor_inplace(grid, fixed);
        int neighborCost = computeCost(grid);

        int delta = neighborCost - currentCost;
        double probability = exp(-delta / T);

        if (delta < 0 || ((double)rand() / RAND_MAX) < probability) {
            // neighbor is better or accepted worse
            currentCost = neighborCost;

            if (currentCost < bestCost) {
                // updating the local best
                bestCost = currentCost;
                for (int i = 0; i < SIZE; i++) {
                    for (int j = 0; j < SIZE; j++) {
                        bestGrid[i][j] = grid[i][j];
                    }
                }
            }
        } 
        else {
            // reverting to the backed-up grid
            for (int i = 0; i < SIZE; i++) {
                for (int j = 0; j < SIZE; j++) {
                    grid[i][j] = backup[i][j];
                }
            }
        }

        T *= alpha; // cooldown
        iteration++;

        if (iteration % 1000 == 0){
            printf("iteration: %d | temp: %.2f | cost: %d\n", iteration, T, bestCost);
        }
        
    }

    // assigning the bestGrid to the original :)
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            grid[i][j] = bestGrid[i][j];
        }
    }

    printf("final cost: %d\n", bestCost);
    printf("iteration: %d\n", iteration);
}