#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "sudoku.h"

GameState gameState;


// dynamic memory allocation since the user can choose the size of the board
int **initializeBoard() {
    int ** newBoard = (int **)malloc(gameState.sizeOfTheBoard * sizeof(int *));
  
    if (newBoard == NULL) {
      perror("FAILED TO ALLOCATE BOARD ROWS");
      exit(EXIT_FAILURE);
    }
  
    for (int i = 0; i < gameState.sizeOfTheBoard; i++) {
      newBoard[i] = (int *)malloc(gameState.sizeOfTheBoard * sizeof(int));
  
      if (newBoard[i] == NULL) {
        perror("FAILED TO ALLOCATE BOARD COLUMNS");
        for (int j = 0; j < i; j++) {
          free(newBoard[j]);
        }
        free(newBoard);
        exit(EXIT_FAILURE);
      }
  
      for (int j = 0; j < gameState.sizeOfTheBoard; j++) {
        newBoard[i][j] = 0;
      }
  
    }
    return newBoard;
  }
  
  
  void freeBoard(){
    if (gameState.board != NULL) {
      for (int i = 0; i < gameState.sizeOfTheBoard; i++) {
        free(gameState.board[i]);
      }
      free(gameState.board);
      gameState.board = NULL; 
    }
  }
  
  
  // printing the board
  void printBoard() {
  
    for (int i = 0; i < gameState.sizeOfTheBoard; i++) {
  
      if (i % gameState.boxSize == 0) {
        printf("--------------------------------\n");
      }
  
      for (int j = 0; j < gameState.sizeOfTheBoard; j++) {
  
        if (j % gameState.boxSize == 0) {
          printf("| ");
        }
  
        if (gameState.board[i][j] == 0) {
          printf(" . ");
        } 
        else if (gameState.board[i][j] <= 9) {
          printf("%2d ", gameState.board[i][j]);
        } 
        else {
          // it will replace 10-16 with a char of A-F
          printf(" %c ", 'A' + (gameState.board[i][j] - 10));
        }
      }
      printf("\n");
    }
  }
  
  
  // SUDOKU FUNCTIONS
  // returns false if 3x3 box contains a number
  bool unUsedInBox(int **grid, const int row, const int col, const int number) {
  
    for (int i = 0; i < gameState.boxSize; i++) {
      for (int j = 0; j < gameState.boxSize; j++) {
        if (grid[row + i][col + j] == number) {
          return false;
        }
      }
    }
    return true;
  }
  
  
  // filling 3x3, random numbers to subgrid
  void fillBox(int **grid, const int row, const int col) {
    int number;
    for (int i = 0; i < gameState.boxSize; i++) {
      for (int j = 0; j < gameState.boxSize; j++) {
        do {
          number = (rand() % gameState.sizeOfTheBoard) + 1;
        } while (!unUsedInBox(grid, row, col, number));
  
        grid[row + i][col + j] = number;
      }
    }
  }
  
  
  // check if it's safe to put a number in row i
  bool unUsedInRow(int **grid, const int i, const int number) {
    for (int j = 0; j < gameState.sizeOfTheBoard; j++) {
      if (grid[i][j] == number) {
        return false;
      }
    }
    return true;
  }
  
  
  // check if it's safe to put a number in column j
  bool unUsedInCol(int **grid, const int j, const int number) {
    for (int i = 0; i < gameState.sizeOfTheBoard; i++) {
      if (grid[i][j] == number) {
        return false;
      }
    }
    return true;
  }
  
  
  // combines all the previous functions to return a bool whether it's safe to put a number in the cell (i, j)
  // ensures that the number isn't used in: row, column, box
  bool checkIfSafe(int **grid, const int i, const int j, const int number) {
    return (unUsedInRow(grid, i, number) && unUsedInCol(grid, j, number) && unUsedInBox(grid, i - i % gameState.boxSize, j - j % gameState.boxSize, number));
  }
  
  
  // filling the diagonal
  void fillDiagonal(int **grid) {
    for (int i = 0; i < gameState.sizeOfTheBoard; i+= gameState.boxSize) {
        fillBox(grid, i, i);
    }
  }
  
  
  // filling the rest with the board
  bool fillRemaining(int **grid, const int i, const int j) {
    // if end of grid
    if (i == gameState.sizeOfTheBoard) {
      return true;
    }
  
    // going to the next row if j == sizeOfTheBoard
    if (j == gameState.sizeOfTheBoard) {
      return fillRemaining(grid, i+1, 0);
    }
  
    //skip if cell is filled (fill diagonal)
    if (grid[i][j] != 0) {
      return fillRemaining(grid, i, j+1);
    }
  
    for (int num = 1; num <= gameState.sizeOfTheBoard; num++) {
      if (checkIfSafe(grid, i, j, num)) {
        grid[i][j] = num;
        if (fillRemaining(grid, i, j+1)) {
          return true;
        }
        grid[i][j] = 0;
      }
    }
    return false;
  }
  
  
  // remove digits randomly
  void removeCells(int **grid) {
  
    int total_cells = gameState.sizeOfTheBoard * gameState.sizeOfTheBoard;
    int cells_to_remove;
  
    switch (gameState.difficulty) {
      case 1:      // EASY
          cells_to_remove = total_cells * 20 / 100;
          break;
      case 2:      // MODERATE
          cells_to_remove = total_cells * 35 / 100;
          break;
      case 3:      // HARD
          cells_to_remove = total_cells * 60 / 100;
          break;
    }
  
    while (cells_to_remove > 0) {
      const int cell = rand() % (gameState.sizeOfTheBoard * gameState.sizeOfTheBoard);
  
      const int i = cell / gameState.sizeOfTheBoard;
      const int j = cell % gameState.sizeOfTheBoard;
  
      if (grid[i][j] != 0) {
        grid[i][j] = 0;
        cells_to_remove--;
      }
    }
  
  }
  
  
  // generating the sudoku board
  int **sudokuGenerator() {
    int ** grid = initializeBoard();
  
    fillDiagonal(grid);
    fillRemaining(grid, 0, 0);
  
    removeCells(grid);
  
    return grid;
  }
  
  
  // empty spots to fill on the board -> spots that the user will have to fill
  int countEmptyCells(int **grid) {
    int count = 0;
    for (int i = 0; i < gameState.sizeOfTheBoard; i++) {
      for (int j = 0; j < gameState.sizeOfTheBoard; j++) {
        if (grid[i][j] == 0) count++;
      }
    }
    return count;
  }
  
  