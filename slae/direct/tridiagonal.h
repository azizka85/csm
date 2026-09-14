#ifndef SLAE_DIRECT_TRIDIAGONAL_H
#define SLAE_DIRECT_TRIDIAGONAL_H

#include <span>

using namespace std;

namespace SLAE::Direct::Tridiagonal {
	struct SolutionState {
		span<double> u;
		span<double> d;
	};

	struct ResidualState {
		span<const double> u;
		span<double> d;
	};

	class ISLAE {
		public:
			virtual void solve(SolutionState state) = 0;
			virtual void residual(ResidualState state) = 0;
	};

	struct UniformMatrixParams {
		double l;
		double l1;

		double c;
		double c0;
		double c1;

		double r;
		double r0;
	};

	class UniformMatrix: ISLAE {
		private:
			UniformMatrixParams params;

		public:
			UniformMatrix(UniformMatrixParams params);

			void solve(SolutionState state) override;
			void residual(ResidualState state) override;

	};

	struct MatrixParams {
		span<const double> l;
		span<const double> c;
		span<const double> r;
	};

	class Matrix : ISLAE {
		private:
			MatrixParams params;

		public:
			Matrix(MatrixParams params);

			void solve(SolutionState state) override;
			void residual(ResidualState state) override;
	};
}

#endif