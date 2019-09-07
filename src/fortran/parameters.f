c     Subroutines for accessing the values from parameters.h, primarily so that
c     they can be used by other languages calling BIRRP subroutines (e.g. to
c     make sure they initialize statically sized arrays to the correct sizes).
c
c     Author: Taylor Viti
      
      subroutine get_npcsm(out)
      integer out
      include "./parameters.h"
      out = npcsm
      return
      end

      subroutine get_nptsm(out)
      integer out
      include "./parameters.h"
      out = nptsm
      return
      end

      subroutine get_nptssm(out)
      integer out
      include "./parameters.h"
      out = nptssm
      return
      end

      subroutine get_noutm(out)
      integer out
      include "./parameters.h"
      out = noutm
      return
      end

      subroutine get_ninpm(out)
      integer out
      include "./parameters.h"
      out = ninpm
      return
      end

      subroutine get_nrefm(out)
      integer out
      include "./parameters.h"
      out = nrefm
      return
      end

      subroutine get_nrsitem(out)
      integer out
      include "./parameters.h"
      out = nrsitem
      return
      end

      subroutine get_nsectm(out)
      integer out
      include "./parameters.h"
      out = nsectm
      return
      end
