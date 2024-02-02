
#include "unit_test.h"

int main(void) {
  Test complex_rigid[] = {
    { .name = "SPC/E Water xy plane, O-H1 on x",
      .test_type = RIGID,
      .atoms = {
        .n = 3,
        .mass = DBLS{15.999, 1.008, 1.008},
        .pos = DBLS{
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          cos(109.47*DEG), sin(109.47*DEG), 0.0,
        },
      },
      .dof = (double[][3]){
        {0.9110346288752766, 0.8995644606470596, 1.0},
        {0.0573987690422075, 0.5373016861966246, 1.0},
        {0.4674075380780282, 0.1272929171608038, 1.0},
      },
    },
  };

  return run_testset((TestSet){
      .name = "Complex Rigid",
      .tests = complex_rigid,
      .n = sizeof(complex_rigid)/sizeof(Test),
    });
}
