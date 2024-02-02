
#include "unit_test.h"

int main(void) {
  Test complex_rigid[] = {
    { .name = "H2",
      .test_type = FLEX,
      .atoms = {
        .n = 2,
        .mass = DBLS{1.008, 1.008},
        .pos = DBLS{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
        },
      },
      .bonds = {
        .n = 1,
        .bonds = (Bond[]){
          {0, 1},
        },
      },
      .dof = (double[][3]){
        {0.5, 1.0, 1.0},
        {0.5, 1.0, 1.0},
      },
    },

    { .name = "Flexible SPC/E Water xy plane, O-H1 on x, 90deg",
      .test_type = FLEX,
      .atoms = {
        .n = 3,
        .mass = DBLS{15.999, 1.008, 1.008},
        .pos = DBLS{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
        },
      },
      .bonds = {
        .n = 2,
        .bonds = (Bond[]){
          {0, 1},
          {0, 2},
        },
      },
      .dof = (double[][3]){
        {0.9407302875286649, 0.9407302875286649, 1.0},
        {0.0592697124713353, 1.0, 1.0},
        {1.0, 0.0592697124713353, 1.0},
      },
    },

    { .name = "Flexible SPC/E Water xy plane, O-H1 on x, 180deg",
      .test_type = FLEX,
      .atoms = {
        .n = 3,
        .mass = DBLS{15.999, 1.008, 1.008},
        .pos = DBLS{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          -1.0, 0.0, 0.0,
        },
      },
      .bonds = {
        .n = 2,
        .bonds = (Bond[]){
          {0, 1},
          {0, 2},
        },
      },
      .dof = (double[][3]){
        {0.8880932556203165, 1.0, 1.0},
        {0.0559533721898418, 1.0, 1.0},
        {0.0559533721898418, 1.0, 1.0},
      },
    },

    { .name = "Flexible Water, 90deg, H as base",
      .test_type = FLEX,
      .atoms = {
        .n = 3,
        .mass = DBLS{1.008, 15.999, 1.008},
        .pos = DBLS{
          1.0, 0.0, 0.0,
          0.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
        },
      },
      .bonds = {
        .n = 2,
        .bonds = (Bond[]){
          {1, 0},
          {1, 2},
        },
      },
      .dof = (double[][3]){
        {0.0592697124713353, 1.0, 1.0},
        {0.9407302875286649, 0.9407302875286649, 1.0},
        {1.0, 0.0592697124713353, 1.0},
      },
    },

    { .name = "Flexible Water, 180deg, sorting needed",
      .test_type = FLEX,
      .atoms = {
        .n = 3,
        .mass = DBLS{15.999, 1.008, 1.008},
        .pos = DBLS{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          -1.0, 0.0, 0.0,
        },
      },
      .bonds = {
        .n = 2,
        .bonds = (Bond[]){
          {0, 2},
          {0, 1},
        },
      },
      .dof = (double[][3]){
        {0.8880932556203165, 1.0, 1.0},
        {0.0559533721898418, 1.0, 1.0},
        {0.0559533721898418, 1.0, 1.0},
      },
    },
  };

  return run_testset((TestSet){
      .name = "Complex Rigid",
      .tests = complex_rigid,
      .n = ntests(complex_rigid),
    });
}
