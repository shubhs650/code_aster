! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine te0542(option, nomte)
!
!
    use resi_refe_module, only: RESI_REFE
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/ltequa.h"
#include "asterfort/nbsigm.h"
#include "asterfort/teattr.h"
#include "asterfort/xbsig.h"
#include "asterfort/xbsir.h"
#include "asterfort/xteddl.h"
#include "asterfort/xteini.h"
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
! FONCTION REALISEE:  CALCUL DES OPTION FORC_NODA ET REFE_FORC_NODA
!                     POUR LES ÉLÉMENTS MECA X-FEM
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfde, jgano, igeom, ivectu
    integer(kind=8) :: jpintt, jcnset, jheavt, jlonch, jbaslo, jvSief, jlsn, jlst
    integer(kind=8) :: jpmilt, ddlm, nfiss, jfisno, jheavn
    integer(kind=8) :: nfh, ddlc, nfe, ibid, ddls, nbsig, nddl, jstno, imate
    integer(kind=8) :: contac, nnom, singu
    aster_logical :: lbid
    real(kind=8) :: sigref(1), depref
    character(len=8) :: enr, elref
    type(RESI_REFE):: refe
! DEB ------------------------------------------------------------------
!
! ---- CARACTERISTIQUES DU TYPE D'ELEMENT :
! ---- GEOMETRIE ET INTEGRATION
!      ------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    call elref1(elref)
!     INITIALISATION DES DIMENSIONS DES DDLS X-FEM
    call xteini(nomte, nfh, nfe, singu, ddlc, &
                nnom, ddls, nddl, ddlm, nfiss, &
                contac)
!
! ---- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!      -----------------------------------------
    nbsig = nbsigm()
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PVECTUR', 'E', ivectu)
!
!     PARAMÈTRES PROPRES À X-FEM
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PHEAVTO', 'L', jheavt)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PBASLOR', 'L', jbaslo)
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
!     PROPRE AUX ELEMENTS 1D ET 2D (QUADRATIQUES)
    call teattr('S', 'XFEM', enr, ibid)
    if (enr(1:2) .eq. 'XH') call jevech('PHEA_NO', 'L', jheavn)
    if ((ibid .eq. 0) .and. ltequa(elref, enr)) &
        call jevech('PPMILTO', 'L', jpmilt)
    if (nfiss .gt. 1) call jevech('PFISNO', 'L', jfisno)
!
    if (option .eq. 'FORC_NODA') then
!      --------------------
! VECTEUR SECOND MEMBRE DONNE EN ENTREE
        call jevech('PSIEFR', 'L', jvSief)
        call jevech('PSTANO', 'L', jstno)
        if (nfe .gt. 0) then
            call jevech('PMATERC', 'L', imate)
        end if
!       CALCUL DU VECTEUR DES FORCES INTERNES (BT*SIGMA)
        call xbsig(ndim, nno, nfh, nfe, ddlc, &
                   ddlm, igeom, jpintt, zi(jcnset), &
                   zi(jheavt), zi(jlonch), zr(jbaslo), zr(jvSief), nbsig, &
                   zr(jlsn), zr(jlst), ivectu, jpmilt, &
                   nfiss, jheavn, jstno, imate)
!
        call xteddl(ndim, nfh, nfe, ddls, nddl, &
                    nno, nnos, zi(jstno), .false._1, lbid, &
                    option, nomte, ddlm, nfiss, jfisno, &
                    vect=zr(ivectu))
!
    else if (option .eq. 'REFE_FORC_NODA') then
!
! --- ON RECUPERE CONTRAINTE ET SAUT DE DEPLACEMENT DE REFERENCE
!
        call refe%Init(nomte)
        sigref(1) = refe%GetRef('SIGM')
        depref = refe%GetRef('DEPL')
        call refe%Check()

        if (nfe .gt. 0) then
            call jevech('PSTANO', 'L', jstno)
        end if
!
! --- ON COMMENCE PAR CALCULER LES CONTRIBUTIONS VOLUMIQUES
!
        call xbsir(ndim, nno, nfh, nfe, ddlc, &
                   ddlm, igeom, jpintt, zi(jcnset), &
                   zi(jheavt), zi(jlonch), zr(jbaslo), sigref, nbsig, &
                   zr(jlsn), zr(jlst), ivectu, jpmilt, &
                   nfiss, jheavn, jstno)
    else
        ASSERT(.false.)
    end if
! FIN ------------------------------------------------------------------
!
end subroutine
