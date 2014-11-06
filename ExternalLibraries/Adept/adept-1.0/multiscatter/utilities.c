#include "ms.h"

/* Print statistics to a stream */
int
ms_print_stats(ms_config config, FILE *file)
{
  fprintf(file, "Single-scattering multiplier: %g\n", config.ss_multiplier);
  fprintf(file, "Total 2-stream source: %g\n", config.total_src);
  fprintf(file, "Total 2-stream reflected: %g\n", config.total_reflected);
  return MS_SUCCESS;
}
