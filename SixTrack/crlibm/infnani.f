      program infnan
************************************************************************
*     (Infinity/NaN)
*     Program to compute and display hexadecimal values of signed NaN
*     and Infinity, and to display hexadecimal and decimal results of
*     their conversion to integer values.
*     (30-Nov-2001)
************************************************************************
      real*4 s
      real*8 d
      real*16 q

      write (6,
     x   '(''Test of signed Infinity and NaN conversion to integers'')')
      write (6,'()')

      s = (1.0 - 1.0)
      write (6,'(a)') 'single precision:'
      write (6,'(8x,a,30x,z8.8)') '+NaN', (1.0 - 1.0) / s
      write (6,'(8x,a,30x,z8.8)') '-NaN', -(1.0 - 1.0) / s
      write (6,'(8x,a,30x,z8.8)') '+Inf', (1.0 + 1.0) / s
      write (6,'(8x,a,30x,z8.8)') '-Inf', -(1.0 + 1.0) / s

      write (6,'(8x,a,17x,z16.8)')
     x    'int(+NaN)', int((1.0 - 1.0) / s)
      write (6,'(8x,a,17x,i16)')
     x    'int(+NaN)', int((1.0 - 1.0) / s)
      write (6,'(8x,a,17x,z16.8)')
     x    'int(-NaN)', int(-(1.0 - 1.0) / s)
      write (6,'(8x,a,17x,i16)')
     x    'int(-NaN)', int(-(1.0 - 1.0) / s)

      write (6,'(8x,a,17x,z16.8)')
     x    'int(+Inf)', int((1.0 + 1.0) / s)
      write (6,'(8x,a,17x,i16)')
     x    'int(+Inf)', int((1.0 + 1.0) / s)
      write (6,'(8x,a,17x,z16.8)')
     x    'int(-Inf)', int(-(1.0 + 1.0) / s)
      write (6,'(8x,a,17x,i16)')
     x    'int(-Inf)', int(-(1.0 + 1.0) / s)

      write (6,'()')

      d = (1.0d+00 - 1.0d+00)
      write (6,'(a)') 'double precision:'

      write (6,'(8x,a,22x,z16.16)')
     x   '+NaN', (1.0d+00 - 1.0d+00) / d
      write (6,'(8x,a,22x,z16.16)')
     x   '-NaN', -(1.0d+00 - 1.0d+00) / d
      write (6,'(8x,a,22x,z16.16)')
     x   '+Inf', (1.0d+00 + 1.0d+00) / d
      write (6,'(8x,a,22x,z16.16)')
     x   '-Inf', -(1.0d+00 + 1.0d+00) / d

      write (6,'(8x,a,17x,z16.8)')
     x    'int(+NaN)', int((1.0d+00 - 1.0d+00) / d)
      write (6,'(8x,a,17x,i16)')
     x    'int(+NaN)', int((1.0d+00 - 1.0d+00) / d)
      write (6,'(8x,a,17x,z16.8)')
     x    'int(-NaN)', int(-(1.0d+00 - 1.0d+00) / d)
      write (6,'(8x,a,17x,i16)')
     x    'int(-NaN)', int(-(1.0d+00 - 1.0d+00) / d)

      write (6,'(8x,a,17x,z16.8)')
     x    'int(+Inf)', int((1.0d+00 + 1.0d+00) / d)
      write (6,'(8x,a,17x,i16)')
     x    'int(+Inf)', int((1.0d+00 + 1.0d+00) / d)
      write (6,'(8x,a,17x,z16.8)')
     x    'int(-Inf)', int(-(1.0d+00 + 1.0d+00) / d)
      write (6,'(8x,a,17x,i16)')
     x    'int(-Inf)', int(-(1.0d+00 + 1.0d+00) / d)
      write (6,'()')

      q = (1.0q+00 - 1.0q+00)
      write (6,'(a)') 'quadruple precision:'

      write (6,'(8x,a,6x,z32.32)')
     x   '+NaN', (1.0q+00 - 1.0q+00) / q
      write (6,'(8x,a,6x,z32.32)')
     x   '-NaN', -(1.0q+00 - 1.0q+00) / q
      write (6,'(8x,a,6x,z32.32)')
     x   '+Inf', (1.0q+00 + 1.0q+00) / q
      write (6,'(8x,a,6x,z32.32)')
     x   '-Inf', -(1.0q+00 + 1.0q+00) / q

      write (6,'(8x,a,17x,z16.8)')
     x    'int(+NaN)', int((1.0q+00 - 1.0q+00) / q)
      write (6,'(8x,a,17x,i16)')
     x    'int(+NaN)', int((1.0q+00 - 1.0q+00) / q)
      write (6,'(8x,a,17x,z16.8)')
     x    'int(+Inf)', int((1.0q+00 + 1.0q+00) / q)
      write (6,'(8x,a,17x,i16)')
     x    'int(+Inf)', int((1.0q+00 + 1.0q+00) / q)
      write (6,'(8x,a,17x,z16.8)')
     x    'int(-NaN)', int((-1.0q+00 - 1.0q+00) / q)
      write (6,'(8x,a,17x,i16)')
     x    'int(-NaN)', int((-1.0q+00 - 1.0q+00) / q)
      write (6,'(8x,a,17x,z16.8)')
     x    'int(-Inf)', int((-1.0q+00 + 1.0q+00) / q)
      write (6,'(8x,a,17x,i16)')
     x    'int(-Inf)', int((-1.0q+00 + 1.0q+00) / q)

      write (6,'()')

      end
