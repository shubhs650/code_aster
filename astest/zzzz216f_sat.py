# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------


from code_aster.Coupling import ExternalCoupling


class FakeSaturne(ExternalCoupling):
    """Fake saturne Solver.

    It allows to test SaturneCoupling without code_saturne installation.

    Arguments:
        debug (bool): Enable debugging mode (default: "False")".
    """

    def __init__(self, debug=False, use_CFEMDEC=False):
        super().__init__("code_saturne", True, debug)
        self._use_CFEMDEC = use_CFEMDEC

    def setup(self, interface, **params):
        """Initialize the coupling.

        Arguments:
            interface (tuple(Mesh, list[str])): whole mesh and groups of the interface.
            params (dict): Parameters of the coupling scheme.
        """

        self.set_parameters(params)

        node_typ = "NODES"
        if self._medcpl.supportOverlap():
            node_typ = "NODES_FE"

        if self._use_CFEMDEC:
            fieldType = "NODES_FE"
        else:
            fieldType = "CELLS"

        self._fields_out = [("boundary_stress", ["FX", "FY", "FZ"], fieldType)]
        self._fields_in = [
            ("mesh_displacement", ["DX", "DY", "DZ"], node_typ),
            ("mesh_velocity", ["DX", "DY", "DZ"], node_typ),
        ]
        self._init_paramedmem(self._other_app, interface)

    def set_parameters(self, params):
        """Set parameters.

        Arguments:
            params (dict): Parameters of the coupling scheme.
        """

        self._params.set_values(params)

        self.MPI.COUPLING_COMM_WORLD.send(0, "ALGOP", int(self._use_CFEMDEC), self.MPI.INT)

        nb_step = (self._params.final_time - self._params.init_time) / self._params.delta_t

        self.MPI.COUPLING_COMM_WORLD.send(0, "NBSSIT", self._params.nb_iter, self.MPI.INT)
        self.MPI.COUPLING_COMM_WORLD.send(0, "TADAPT", int(self._params.adapt_step), self.MPI.INT)

        self.MPI.COUPLING_COMM_WORLD.send(0, "EPSILO", self._params.epsilon, self.MPI.DOUBLE)

        self.MPI.COUPLING_COMM_WORLD.send(0, "TTINIT", self._params.init_time, self.MPI.DOUBLE)
        self.MPI.COUPLING_COMM_WORLD.send(0, "TTMAX", self._params.final_time, self.MPI.DOUBLE)
        self.MPI.COUPLING_COMM_WORLD.send(0, "PDTREF", self._params.delta_t, self.MPI.DOUBLE)

    def run(self, solver):
        """Execute the coupling loop.

        Arguments:
            solver (object): Solver contains at least a method run_iteration.
        """

        stepper = self._params.stepper
        first_start = self._starter
        completed = False
        input_data = None
        istep = 0

        if not self._params.adapt_step:
            self.sync(end_coupling=False)

        while not stepper.isFinished():
            istep += 1
            current_time = stepper.getCurrent()
            delta_time = stepper.getIncrement()

            _ = self.MPI.COUPLING_COMM_WORLD.recv(istep, "DTAST", self.MPI.DOUBLE)
            self.MPI.COUPLING_COMM_WORLD.send(istep, "DTCALC", delta_time, self.MPI.DOUBLE)

            print("coupling step #{0:d}, time: {1:f}".format(istep, current_time))

            for i_iter in range(self._params.nb_iter):

                print("coupling iteration #{0:d}, time: {1:f}".format(i_iter, current_time))

                if first_start:
                    first_start = False
                    input_data = {}
                    for name, _, _ in self._fields_in:
                        input_data[name] = None
                elif i_iter > 0:
                    input_data = self.recv_input_fields()

                has_cvg = solver.run_iteration(i_iter, current_time, delta_time, input_data)

                # send cvg
                self.MPI.COUPLING_COMM_WORLD.send(istep, "ICVAST", int(has_cvg), self.MPI.INT)

                # send results to code_aster
                self.send_output_fields()

                if has_cvg:
                    break

            stepper.completed()
            input_data = self.recv_input_fields()

            if not self._params.adapt_step:
                exit_coupling = self.sync(end_coupling=stepper.isFinished())
            else:
                exit_coupling = stepper.isFinished()

            print(f"end of time step {current_time} with status: {exit_coupling}")

        self.sync(end_coupling=stepper.isFinished())

        print(
            "coupling {0} with exit status: {1}".format(
                "completed" if completed else "interrupted", exit_coupling
            )
        )
