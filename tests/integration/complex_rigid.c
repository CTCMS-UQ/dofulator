#include "integration_test.h"

int main(void) {
  Test complex_rigid[] = {
    { .name = "SPC/E Water xy plane, O-H1 on x",
      .test_type = RIGID,
      .atoms = {
        .n = 3,
        .mass = DBLS{15.999, 1.008, 1.008},
        .x = DBLSX{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          cos(109.47*DEG), sin(109.47*DEG), 0.0,
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
        {0.9110346288752766, 0.8995644606470596, 1.0},
        {0.0573987690422075, 0.5373016861966246, 1.0},
        {0.4674075380780282, 0.1272929171608038, 1.0},
      },
    },
    {
      .name = "Large fragment, merging required",
      .test_type = RIGID,
      .atoms = {
        .n = 5,
        .mass = DBLS{1., 1., 1., 1., 1.},
        .x = DBLSX{
          0., 0., 0.,
          1., 0., 0.,
          2., 0., 0.,
          3., 0., 0.,
          4., 0., 0.,
        },
      },
      .bonds = {
        .n = 4,
        .bonds = (Bond[]){
          {0, 1},
          {2, 3},   // Create 2nd fragment to be merged
          {2, 1},   // Merge two fragments
          {3, 4},   // Find frag 0 through invalidated frag 1
        },
      },
      .dof = (double[][3]){
        {0.2, 0.6, 0.6},
        {0.2, 0.3, 0.3},
        {0.2, 0.2, 0.2},
        {0.2, 0.3, 0.3},
        {0.2, 0.6, 0.6},
      },
    },
    { .name = "Multiple fragments",
      .test_type = RIGID,
      .atoms = {
        .n = 8,
        .mass = DBLS{15.999, 1.008, 1.008, 1., 1., 1., 1., 1.},
        .x = DBLSX{
          // H2O
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          cos(109.47*DEG), sin(109.47*DEG), 0.0,
          // Dumbbell
          5., 1., 2.,
          6., 1., 2.,
          // Linear 3 atom
          1., -5., -3.,
          1., -4., -3.,
          1., -3., -3.,
        },
      },
      .bonds = {
        .n = 5,
        .bonds = (Bond[]){
          {0, 1},
          {0, 2},
          {3, 4},
          {6, 7},
          {5, 6},
        },
      },
      .dof = (double[][3]){
        // H2O
        {0.9110346288752766, 0.8995644606470596, 1.0},
        {0.0573987690422075, 0.5373016861966246, 1.0},
        {0.4674075380780282, 0.1272929171608038, 1.0},
        // Dumbbell
        {0.5, 1.0, 1.0},
        {0.5, 1.0, 1.0},
        // 3 atom linear
        {0.8333333333333333, 0.3333333333333333, 0.8333333333333333},
        {0.3333333333333333, 0.3333333333333333, 0.3333333333333333},
        {0.8333333333333333, 0.3333333333333333, 0.8333333333333333},
      },
    },
  };

  return run_testset((TestSet){
      .name = "Complex Rigid",
      .tests = complex_rigid,
      .n = ntests(complex_rigid),
    });
}
