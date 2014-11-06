/* config.c -- Configuration manipulation for multiscatter algorithms

   Copyright (C) 2007 Robin Hogan <r.j.hogan@reading.ac.uk> 

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "ms.h"

void
ms_set_default_config(ms_config* config)
{
  static const ms_config local_config = MS_DEFAULT_CONFIG;
  *config = local_config;
}

void 
ms_set_options(ms_config* config, int new_options)
{
  config->options = new_options;
}
void
ms_add_options(ms_config* config, int new_options)
{
  config->options |= new_options;
}


int
ms_set_max_scattering_order(ms_config* config, int norder)
{
  if (norder > 1) {
    config->max_scattering_order = norder;
    return MS_SUCCESS;
  }
  else {
    return MS_FAILURE;
  }
}


void
ms_set_first_wide_angle_gate(ms_config* config, int first_gate)
{
  if (config) {
    config->first_wide_angle_gate = first_gate;
  }
}
