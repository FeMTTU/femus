// simulate_radiances.h - a function taking inactive arguments that returns also Jacobian matrices

void simulate_radiances(int n, // Size of temperature array
			// Input variables:
			double surface_temperature, 
			const double* temperature,
			// Output variables:
			double radiance[2],
			// Output Jacobians:
			double dradiance_dsurface_temperature[2],
			double* dradiance_dtemperature);
