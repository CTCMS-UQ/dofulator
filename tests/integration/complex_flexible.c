#include "integration_test.h"

int main(void) {
  Test basic_flexible[] = {
    { .name = "Square (single loop closure)",
      .test_type = FLEX,
      .atoms = {
        .n = 4,
        .mass = DBLS{1., 1., 1., 1.},
        .x = DBLSX{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          1.0, 1.0, 0.0,
          0.0, 1.0, 0.0,
        },
      },
      .bonds = {
        .n = 4,
        .bonds = (Bond[]){
          {0, 1},
          {1, 2},
          {2, 3},
          {3, 0},
        },
      },
      .dof = (double[][3]){
        {0.5, 0.5, 1.0},
        {0.5, 0.5, 1.0},
        {0.5, 0.5, 1.0},
        {0.5, 0.5, 1.0},
      },
    },

    { .name = "Cube (multiple loop closures)",
      .test_type = FLEX,
      .atoms = {
        .n = 8,
        .mass = DBLS{1., 1., 1., 1., 1., 1., 1., 1.},
        .x = DBLSX{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          1.0, 1.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          1.0, 0.0, 1.0,
          1.0, 1.0, 1.0,
          0.0, 1.0, 1.0,
        },
      },
      .bonds = {
        .n = 12,
        .bonds = (Bond[]){
          // Bottom square
          {0, 1},
          {1, 2},
          {2, 3},
          {3, 0},
          // Sides
          {0, 4},
          {1, 5},
          {2, 6},
          {3, 7},
          // Top square
          {4, 5},
          {5, 6},
          {6, 7},
          {7, 4},
        },
      },
      .dof = (double[][3]){
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
      },
    },

    { .name = "Cube (loop closures + merge with loops)",
      .test_type = FLEX,
      .atoms = {
        .n = 8,
        .mass = DBLS{1., 1., 1., 1., 1., 1., 1., 1.},
        .x = DBLSX{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          1.0, 1.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          1.0, 0.0, 1.0,
          1.0, 1.0, 1.0,
          0.0, 1.0, 1.0,
        },
      },
      .bonds = {
        .n = 12,
        .bonds = (Bond[]){
          // Bottom square
          {0, 1},
          {1, 2},
          {2, 3},
          {3, 0},
          // Top square
          {4, 5},
          {5, 6},
          {6, 7},
          {7, 4},
          // Sides
          {0, 4},
          {1, 5},
          {2, 6},
          {3, 7},
        },
      },
      .dof = (double[][3]){
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5},
      },
    },
  };

  return run_testset((TestSet){
      .name = "Complex Flexible",
      .tests = basic_flexible,
      .n = ntests(basic_flexible),
    });
}
