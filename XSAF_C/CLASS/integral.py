from scipy.interpolate import CubicSpline
import numpy as np


def get_integral(pk, kh):
    P_m = CubicSpline(np.log10(kh[-1]), np.log10(pk[-1]))
    integrale = (
        1.0
        / (2 * np.pi ** 2)
        * np.sum(
            (kh[-1][1:] - kh[-1][:-1])
            / 6.0
            * (
                (
                    kh[-1][:-1] ** 2
                    * (
                        3
                        * (
                            np.sin(8 * kh[-1][:-1])
                            - 8 * kh[-1][:-1] * np.cos(8 * kh[-1][:-1])
                        )
                        / (8 * kh[-1][:-1]) ** 3
                    )
                    ** 2
                    * pk[-1][:-1]
                )
                + (
                    kh[-1][1:] ** 2
                    * (
                        3
                        * (
                            np.sin(8 * kh[-1][1:])
                            - 8 * kh[-1][1:] * np.cos(8 * kh[-1][1:])
                        )
                        / (8 * kh[-1][1:]) ** 3
                    )
                    ** 2
                    * pk[-1][1:]
                )
                + 4.0
                * (
                    ((kh[-1][1:] + kh[-1][:-1]) / 2.0) ** 2
                    * (
                        3
                        * (
                            np.sin(8 * (kh[-1][1:] + kh[-1][:-1]) / 2.0)
                            - 8
                            * (kh[-1][1:] + kh[-1][:-1])
                            / 2
                            * np.cos(8 * (kh[-1][1:] + kh[-1][:-1]) / 2.0)
                        )
                        / (8 * (kh[-1][1:] + kh[-1][:-1]) / 2.0) ** 3
                    )
                    ** 2
                    * 10 ** P_m(np.log10((kh[-1][:-1] + kh[-1][1:]) / 2))
                )
            )
        )
    )
    return integrale
