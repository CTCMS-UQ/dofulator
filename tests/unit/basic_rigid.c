#include "unit_test.h"

int main(void) {
  Test basic_rigid[] = {
    { .name = "Single atom, origin",
      .test_type = RIGID,
      .atoms = {
        .n = 1,
        .mass = DBLS{1.0},
        .pos =  DBLS{0.0, 0.0, 0.0},
      },
      .dof = (double[][3]){
        {1.0, 1.0, 1.0},
      },
    },

    { .name = "Single atom, translated",
      .test_type = RIGID,
      .atoms = {
        .n = 1,
        .mass = DBLS{1.0},
        .pos =  DBLS{1.0, 3.0, -2.0},
      },
      .dof = (double[][3]){
        {1.0, 1.0, 1.0},
      },
    },

    { .name = "Dumbbell x\0",
      .test_type = RIGID,
      .atoms = {
        .n = 2,
        .mass = DBLS{1.0, 1.0},
        .pos =  DBLS{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0
        },
      },
      .dof = (double[][3]){
        {0.5, 1.0, 1.0},
        {0.5, 1.0, 1.0},
      },
    },

    { .name = "Dumbbell y",
      .test_type = RIGID,
      .atoms = {
        .n = 2,
        .mass = DBLS{1.0, 1.0},
        .pos =  DBLS{
          0.0, 0.0, 0.0,
          0.0, 1.0, 0.0
        },
      },
      .dof = (double[][3]){
        {1.0, 0.5, 1.0},
        {1.0, 0.5, 1.0},
      },
    },

    { .name = "Dumbbell z",
      .test_type = RIGID,
      .atoms = {
        .n = 2,
        .mass = DBLS{1.0, 1.0},
        .pos =  DBLS{
          0.0, 0.0, 0.0,
          0.0, 0.0, 1.0
        },
      },
      .dof = (double[][3]){
        {1.0, 1.0, 0.5},
        {1.0, 1.0, 0.5},
      },
    },

    { .name = "Cube",
      .test_type = RIGID,
      .atoms = {
        .n = 8,
        .mass = DBLS{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
        .pos =  DBLS{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          1.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          1.0, 0.0, 1.0,
          0.0, 1.0, 1.0,
          1.0, 1.0, 1.0,
        },
      },
      .dof = (double[][3]){
        {0.25, 0.25, 0.25},
        {0.25, 0.25, 0.25},
        {0.25, 0.25, 0.25},
        {0.25, 0.25, 0.25},
        {0.25, 0.25, 0.25},
        {0.25, 0.25, 0.25},
        {0.25, 0.25, 0.25},
        {0.25, 0.25, 0.25},
      }
    },
  };

  return run_testset((TestSet){
      .name = "Basic Rigid",
      .tests = basic_rigid,
      .n = sizeof(basic_rigid) / sizeof(Test),
    });
}
