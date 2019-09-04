!
!  Copyright 2017-2019 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!

subroutine hamiltonian(zu,flag_current)
  use Global_Variables
  implicit none
  complex(8), intent(inout) :: zu(NL,NBoccmax,NK_s:NK_e)
  logical, intent(in)       :: flag_current
  integer, parameter :: ID_PROPAGATOR_TAYLOR = 1
  integer, parameter :: ID_PROPAGATOR_LANCZOS = 2
  integer, parameter :: ID_PROPAGATOR = ID_PROPAGATOR_TAYLOR

  select case(ID_PROPAGATOR)
  case(ID_PROPAGATOR_TAYLOR)
    call hamiltonian_Taylor(zu,flag_current)
  case(ID_PROPAGATOR_LANCZOS)
    call hamiltonian_Lanczos(zu,flag_current)
  case default
    call err_finalize('invalid propagator')
  end select


end subroutine hamiltonian
