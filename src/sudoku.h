#ifndef SUDOKU_H
#define SUDOKU_H

#include <stdbool.h>

typedef struct {
    int sizeOfTheBoard;
    int boxSize;
    int **board;
    int hearts;
    int correctMoves;
    int difficulty;
} GameState;

extern GameState gameState;

int **initializeBoard();
void freeBoard();
void printBoard();
bool unUsedInBox(int **grid, int row, int col, int number);
void fillBox(int **grid, int row, int col);
bool unUsedInRow(int **grid, int i, int number);
bool unUsedInCol(int **grid, int j, int number);
bool checkIfSafe(int **grid, int i, int j, int number);
void fillDiagonal(int **grid);
bool fillRemaining(int **grid, int i, int j);
void removeCells(int **grid);
int **sudokuGenerator();
int countEmptyCells(int **grid);

#endif
