

      module hydfinal
      !DEC$ ATTRIbUTES DLLEXPORT :: SP, DP, QP, PR
      
      public
   
      integer, parameter              :: SP = KIND(1.0)
      integer, parameter              :: DP = KIND(1.d0)
     ! integer, parameter              :: QP = KIND(1.q0) error 658 - Invalid character '0' in operator name, only 'a-z' and 'A-Z' may be used
      integer, parameter              :: PR = DP
      !real(PR), parameter, private    :: pi = ACOS(-1.0_PR)
      !complex(PR), parameter, private :: pi_c = (pi, 0.0_PR)

      !--- DEFINITION OF THE DERIVED TYPE
      TYPE hyperdual
        double precision :: f0
        double precision :: f1
        double precision :: f2
        double precision :: f12
      END TYPE

      
      !================================================================!
      !                Overloading hyperdual functions                 !
      !================================================================!
      
      !----- Constructor
      interface hyd
        module procedure hyd_from_dble
      end interface

      !----- Assignment
      interface assignment(=)
        module procedure hyd_assign_hyd, hyd_assign_cplx, hyd_assign_dble, hyd_assign_int, &
            hyd_array_assign_hyd, hyd_array_assign_hyd_array, hyd_array_assign_dble, &
              hyd_array_assign_dble_array, hyd_array_assign_int, hyd_array_assign_int_array, &
            hyd_matrix_assign_hyd, hyd_matrix_assign_hyd_matrix, hyd_matrix_assign_dble, &
              hyd_matrix_assign_dble_matrix, hyd_matrix_assign_int, hyd_matrix_assign_int_matrix, &
            hyd_tens_assign_hyd, hyd_tens_assign_hyd_tens, hyd_tens_assign_dble, &
              hyd_tens_assign_dble_tens, hyd_tens_assign_int, hyd_tens_assign_int_tens
      end interface

      !----- Comparison operators      
      interface operator(==)
        module procedure hyd_eq_hyd, hyd_eq_dble, dble_eq_hyd
      end interface

      interface operator(/=)
        module procedure hyd_ne_hyd, hyd_ne_dble, dble_ne_hyd
      end interface

      interface operator(>)
        module procedure hyd_gt_hyd, hyd_gt_dble, dble_gt_hyd
      end interface
      
      interface operator(>=)
        module procedure hyd_ge_hyd, hyd_ge_dble, dble_ge_hyd
      end interface
      
      interface operator(<)
        module procedure hyd_lt_hyd, hyd_lt_dble, dble_lt_hyd
      end interface
      
      interface operator(<=)
        module procedure hyd_le_hyd, hyd_le_dble, dble_le_hyd
      end interface
      

      !----- Arithmetic operators
      interface operator (+)
        module procedure hyd_plus_hyd, hyd_plus_hyd_array, hyd_plus_hyd_matrix, &
              hyd_plus_hyd_tens, hyd_plus_dble, hyd_plus_dble_array, hyd_plus_dble_matrix, &
              hyd_plus_dble_tens, hyd_plus_int, hyd_plus_int_array, &
              hyd_plus_int_matrix, hyd_plus_int_tens, &
            hyd_array_plus_hyd, hyd_array_plus_hyd_array, hyd_array_plus_dble, &
              hyd_array_plus_dble_array, hyd_array_plus_int, &
              hyd_array_plus_int_array, &
            hyd_matrix_plus_hyd, hyd_matrix_plus_hyd_matrix, hyd_matrix_plus_dble, &
              hyd_matrix_plus_dble_matrix, hyd_matrix_plus_int, &
              hyd_matrix_plus_int_matrix, &
            hyd_tens_plus_hyd, hyd_tens_plus_hyd_tens, hyd_tens_plus_dble, &
              hyd_tens_plus_dble_tens, hyd_tens_plus_int, &
              hyd_tens_plus_int_tens, &
            dble_plus_hyd, dble_plus_hyd_array, dble_plus_hyd_matrix, dble_plus_hyd_tens, &
            dble_array_plus_hyd, dble_array_plus_hyd_array, &
            dble_matrix_plus_hyd, dble_matrix_plus_hyd_matrix, &
            dble_tens_plus_hyd, dble_tens_plus_hyd_tens, &
            int_plus_hyd, int_plus_hyd_array, int_plus_hyd_matrix, &
              int_plus_hyd_tens, &
            int_array_plus_hyd, int_array_plus_hyd_array, &
            int_matrix_plus_hyd, int_matrix_plus_hyd_matrix, &
            int_tens_plus_hyd, int_tens_plus_hyd_tens
      end interface

      interface operator(-)
         module procedure hyd_minus_hyd, hyd_minus_hyd_array, hyd_minus_hyd_matrix, &
             hyd_minus_hyd_tens, hyd_minus_dble, hyd_minus_dble_array, hyd_minus_dble_matrix,&
             hyd_minus_dble_tens, hyd_minus_int, hyd_minus_int_array, hyd_minus_int_matrix, &
             hyd_minus_int_tens, &
            hyd_array_minus_hyd, hyd_array_minus_hyd_array, &
              hyd_array_minus_dble, hyd_array_minus_dble_array, &
              hyd_array_minus_int, hyd_array_minus_int_array, &
            hyd_matrix_minus_hyd, hyd_matrix_minus_hyd_matrix, &
              hyd_matrix_minus_dble, hyd_matrix_minus_dble_matrix, &
              hyd_matrix_minus_int, hyd_matrix_minus_int_matrix, &
            hyd_tens_minus_hyd, hyd_tens_minus_hyd_tens, &
              hyd_tens_minus_dble, hyd_tens_minus_dble_tens, &
              hyd_tens_minus_int, hyd_tens_minus_int_tens, &
            dble_minus_hyd, dble_minus_hyd_array, &
              dble_minus_hyd_matrix, dble_minus_hyd_tens, &
            dble_array_minus_hyd, dble_array_minus_hyd_array, &
            dble_matrix_minus_hyd, dble_matrix_minus_hyd_matrix, &
            dble_tens_minus_hyd, dble_tens_minus_hyd_tens, &
            int_minus_hyd, int_minus_hyd_array, &
              int_minus_hyd_matrix, int_minus_hyd_tens, &
            int_array_minus_hyd, int_array_minus_hyd_array, &
            int_matrix_minus_hyd, int_matrix_minus_hyd_matrix, &
            int_tens_minus_hyd, int_tens_minus_hyd_tens, &
            minus_hyd, minus_hyd_array, minus_hyd_matrix, &
              minus_hyd_tens
      end interface

      interface operator(*)
        module procedure hyd_mul_hyd, hyd_mul_hyd_array, hyd_mul_hyd_matrix, &
              hyd_mul_hyd_tens, hyd_mul_dble, hyd_mul_dble_array, hyd_mul_dble_matrix, &
              hyd_mul_dble_tens, hyd_mul_int, hyd_mul_int_array, hyd_mul_int_matrix, &
              hyd_mul_int_tens, & 
            hyd_array_mul_hyd, hyd_array_mul_dble, hyd_array_mul_int, &
            hyd_matrix_mul_hyd, hyd_matrix_mul_dble, hyd_matrix_mul_int, &
            hyd_tens_mul_hyd, hyd_tens_mul_dble, hyd_tens_mul_int, &
            dble_mul_hyd, dble_mul_hyd_array, dble_mul_hyd_matrix, &
              dble_mul_hyd_tens, &
            dble_array_mul_hyd, &
            dble_matrix_mul_hyd, &
            dble_tens_mul_hyd, &
            int_mul_hyd, int_mul_hyd_array, int_mul_hyd_matrix, &
              int_mul_hyd_tens, &
            int_array_mul_hyd, &
            int_matrix_mul_hyd, &
            int_tens_mul_hyd
      end interface

      interface operator(/)
        module procedure hyd_div_hyd, hyd_div_dble, hyd_div_int, &
            hyd_array_div_hyd, hyd_array_div_dble, hyd_array_div_int, &
            hyd_matrix_div_hyd, hyd_matrix_div_dble, hyd_matrix_div_int, &
            hyd_tens_div_hyd, hyd_tens_div_dble, hyd_tens_div_int, &
            dble_div_hyd, int_div_hyd, &
            dble_array_div_hyd, int_array_div_hyd, &
            dble_matrix_div_hyd, int_matrix_div_hyd, &
            dble_tens_div_hyd, int_tens_div_hyd
      end interface

      interface operator(**)
        module procedure hyd_pow_hyd, hyd_pow_int, hyd_pow_dble
      end interface

      !----- Intrinsic functions
      interface dble
        module procedure dble_hyd, dble_hyd_array, dble_hyd_matrix
      end interface
      
      interface abs
        module procedure abs_hyd, abs_hyd_array, abs_hyd_matrix
      end interface
      
      interface sign
        module procedure sign_hyd_hyd
      end interface

      interface max
        module procedure max_hyd_hyd, max_hyd_dble, max_dble_hyd
      end interface

      interface min
        module procedure min_hyd_hyd, min_hyd_dble, min_dble_hyd
      end interface

      interface maxval
        module procedure maxval_hyd_array
      end interface

      interface cos
        module procedure cos_hyd
      end interface
      
      interface sin
        module procedure sin_hyd
      end interface

      interface tan
        module procedure tan_hyd
      end interface
      
      interface sqrt
        module procedure sqrt_hyd
      end interface

      interface acos
        module procedure acos_hyd
      end interface

      interface asin
        module procedure asin_hyd
      end interface

      interface atan
        module procedure atan_hyd
      end interface
      
      interface matmul
        module procedure &
          matmul_hyd_array_hyd_matrix, matmul_hyd_array_dble_matrix, &
          matmul_hyd_matrix_hyd_array, matmul_hyd_matrix_hyd_matrix, &
            matmul_hyd_matrix_dble_array, matmul_hyd_matrix_dble_matrix, &
          matmul_dble_array_hyd_matrix, &
          matmul_dble_matrix_hyd_array, matmul_dble_matrix_hyd_matrix
      end interface

      interface dot_product
        module procedure &
          dot_product_hyd_array_hyd_array, dot_product_hyd_array_dble_array, &
          dot_product_dble_array_hyd_array
      end interface
      
      interface exp
        module procedure exp_hyd
      end interface

      interface log
        module procedure log_hyd
      end interface
      
      interface aimag
          module procedure imag1_hyd, imag2_array_hyd
      end interface

      interface conjg
          module procedure conjg_hyd
      end interface

      
      !================================================================!
      !                 Overloading COMPLEX functions                  !
      !================================================================!
 



    

      
    
!==========================================================================!
!==========================================================================!

      CONTAINS

      !================================================================!
      !                Overloading hyperdual functions                 !                                                                               h
      !================================================================!

      !----------------------------------------------------------------!
      !                         CONSTRUCTOR                            !
      !----------------------------------------------------------------!

        
        function hyd_from_dble(f0, f1, f2, f12) result(q)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_from_dble" :: hyd_from_dble
                
          double precision, intent(in) :: f0, f1, f2, f12
          TYPE(hyperdual)         :: q 
                
          q%f0 = f0
          q%f1 = f1
          q%f2 = f2
          q%f12 = f12
                    
          return
        end function hyd_from_dble

      !----------------------------------------------------------------!
      !                          ASSIGNMENT                            !
      !----------------------------------------------------------------!

        subroutine hyd_assign_hyd(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_assign_hyd" :: hyd_assign_hyd

          
          TYPE(hyperdual), intent(out)  :: lhs
          TYPE(hyperdual), intent(in)   :: rhs
        
          lhs%f0 = rhs%f0
          lhs%f1 = rhs%f1
          lhs%f1 = rhs%f2
          lhs%f1 = rhs%f12
        
        end subroutine hyd_assign_hyd

        subroutine hyd_assign_cplx(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_assign_cplx" :: hyd_assign_cplx

          
          TYPE(hyperdual), intent(out)  :: lhs
          complex(PR), intent(in)       :: rhs
        
          lhs%f0 = rhs
          lhs%f1 = CMPLX(0.0_PR, 0.0_PR, PR)
        
        end subroutine hyd_assign_cplx

        subroutine hyd_assign_dble(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_assign_dble" :: hyd_assign_dble
        
                    
          TYPE(hyperdual), intent(out)  :: lhs
          real(PR), intent(in)          :: rhs
        
          lhs%f0 = rhs
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
        
        end subroutine hyd_assign_dble
        
        subroutine hyd_assign_int(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_assign_int" :: hyd_assign_int
        
          
          TYPE(hyperdual), intent(out)  :: lhs
          integer, intent(in)           :: rhs
        
          lhs%f0 = REAL(rhs,PR)
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
        
        end subroutine hyd_assign_int

        subroutine hyd_array_assign_hyd(lhs, rhs)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_assign_hyd" :: hyd_array_assign_hyd

          
          TYPE(hyperdual), dimension(:), intent(out)  :: lhs
          TYPE(hyperdual), intent(in)                 :: rhs

          lhs%f0 = rhs%f0
          lhs%f1 = rhs%f1
          lhs%f2 = rhs%f2
          lhs%f12 = rhs%f12

        end subroutine hyd_array_assign_hyd

        subroutine hyd_array_assign_hyd_array(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_assign_hyd_array" :: hyd_array_assign_hyd_array
        
          
          TYPE(hyperdual), dimension(:), intent(out)  :: lhs
          TYPE(hyperdual), dimension(:), intent(in)   :: rhs
                  
          lhs%f0 = rhs%f0
          lhs%f1 = rhs%f1
          lhs%f2 = rhs%f2
          lhs%f12 = rhs%f12
          
        end subroutine hyd_array_assign_hyd_array

        subroutine hyd_array_assign_dble(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_assign_dble" :: hyd_array_assign_dble
                
          
          TYPE(hyperdual), dimension(:), intent(out)  :: lhs
          real(PR), intent(in)                        :: rhs
          lhs%f0 = rhs
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
          
                
        end subroutine hyd_array_assign_dble

        subroutine hyd_array_assign_dble_array(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_assign_dble_array" :: hyd_array_assign_dble_array
            
          
          TYPE(hyperdual), dimension(:), intent(out)  :: lhs
          real(PR), dimension(:), intent(in)          :: rhs
                
          lhs%f0 = rhs
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
                
        end subroutine hyd_array_assign_dble_array

        subroutine hyd_array_assign_int(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_assign_int" :: hyd_array_assign_int
                
          
          TYPE(hyperdual), dimension(:), intent(out)  :: lhs
          integer, intent(in)                         :: rhs
                  
          lhs%f0 = rhs
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
            
        end subroutine hyd_array_assign_int

        subroutine hyd_array_assign_int_array(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_assign_int_array" :: hyd_array_assign_int_array
                
          
          TYPE(hyperdual), dimension(:), intent(out)  :: lhs
          integer, dimension(:), intent(in)           :: rhs
                  
          lhs%f0 = Real(rhs, PR)
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
            
        end subroutine hyd_array_assign_int_array

        subroutine hyd_matrix_assign_hyd(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_assign_hyd" :: hyd_matrix_assign_hyd
        
          
          TYPE(hyperdual), dimension(:,:), intent(out)  :: lhs
          TYPE(hyperdual), intent(in)                   :: rhs
        
          lhs%f0 = rhs%f0
          lhs%f1 = rhs%f1
          lhs%f2 = rhs%f2
          lhs%f12 = rhs%f12
        
        end subroutine hyd_matrix_assign_hyd

        subroutine hyd_matrix_assign_hyd_matrix(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_assign_hyd_matrix" :: hyd_matrix_assign_hyd_matrix
        
          
          TYPE(hyperdual), dimension(:,:), intent(out)  :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in)   :: rhs
        
          lhs%f0 = rhs%f0
          lhs%f1 = rhs%f1
          lhs%f2 = rhs%f2
          lhs%f12 = rhs%f12
        
        end subroutine hyd_matrix_assign_hyd_matrix

        subroutine hyd_matrix_assign_dble(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_assign_dble" :: hyd_matrix_assign_dble
        
          
          TYPE(hyperdual), dimension(:,:), intent(out)  :: lhs
          real(PR), intent(in)                          :: rhs
        
          lhs%f0 = rhs
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
        
        end subroutine hyd_matrix_assign_dble

        subroutine hyd_matrix_assign_dble_matrix(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_assign_dble_matrix" :: hyd_matrix_assign_dble_matrix
                
          
          TYPE(hyperdual), dimension(:,:), intent(out)  :: lhs
          real(PR), dimension(:,:), intent(in)          :: rhs
                
          lhs%f0 = rhs
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
              
        end subroutine hyd_matrix_assign_dble_matrix

        subroutine hyd_matrix_assign_int(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_assign_int" :: hyd_matrix_assign_int
        
          
          TYPE(hyperdual), dimension(:,:), intent(out)  :: lhs
          integer, intent(in)                           :: rhs
        
          lhs%f0 = Real(rhs, PR)
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
          
        end subroutine hyd_matrix_assign_int

        subroutine hyd_matrix_assign_int_matrix(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_assign_int_matrix" :: hyd_matrix_assign_int_matrix
        
          
          TYPE(hyperdual), dimension(:,:), intent(out)  :: lhs
          integer, dimension(:,:), intent(in)           :: rhs
        
          lhs%f0 = Real(rhs, PR)
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
        
        end subroutine hyd_matrix_assign_int_matrix

        subroutine hyd_tens_assign_hyd(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_assign_hyd" :: hyd_tens_assign_hyd
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(out)  :: lhs
          TYPE(hyperdual), intent(in)                     :: rhs
        
          lhs%f0 = rhs%f0
          lhs%f1 = rhs%f1
          lhs%f2 = rhs%f2
          lhs%f12 = rhs%f12
        
        end subroutine hyd_tens_assign_hyd

        subroutine hyd_tens_assign_hyd_tens(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_assign_hyd_tens" :: hyd_tens_assign_hyd_tens
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(out)  :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in)   :: rhs
        
          lhs%f0 = rhs%f0
          lhs%f1 = rhs%f1
          lhs%f2 = rhs%f2
          lhs%f12 = rhs%f12
        
        end subroutine hyd_tens_assign_hyd_tens

        subroutine hyd_tens_assign_dble(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_assign_dble" :: hyd_tens_assign_dble
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(out) :: lhs
          real(PR), intent(in) :: rhs
        
          lhs%f0 = rhs
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
        
        end subroutine hyd_tens_assign_dble

        subroutine hyd_tens_assign_dble_tens(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_assign_dble_tens" :: hyd_tens_assign_dble_tens
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(out)  :: lhs
          real(PR), dimension(:,:,:), intent(in)          :: rhs
        
          lhs%f0 = rhs
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
        
        end subroutine hyd_tens_assign_dble_tens

        subroutine hyd_tens_assign_int(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_assign_int" :: hyd_tens_assign_int
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(out)  :: lhs
          integer, intent(in)                             :: rhs
        
          lhs%f0 = Real(rhs, PR)
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
          
        end subroutine hyd_tens_assign_int

        subroutine hyd_tens_assign_int_tens(lhs, rhs) 
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_assign_int_tens" :: hyd_tens_assign_int_tens
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(out)  :: lhs
          integer, dimension(:,:,:), intent(in)           :: rhs
        
          lhs%f0 = Real(rhs, PR)
          lhs%f1 = 0.0_PR
          lhs%f2 = 0.0_PR
          lhs%f12 = 0.0_PR
          
        end subroutine hyd_tens_assign_int_tens

      !----------------------------------------------------------------!
      !                     COMPARISON OPERATORS                       !
      !----------------------------------------------------------------!
      ! As mentioned in Lantoine [1] and Martins [2], the comparison operators
      ! must be written so that they ensure that the same branch is followed
      ! whether the input is real, complex or hyperdual. They will therefore
      ! only use the real parts of the input: real(cplx) or real(hyd%f0)

      !----- EQ operator (==)
        function hyd_eq_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_eq_hyd" :: hyd_eq_hyd

          
          TYPE(hyperdual), intent(in) :: lhs, rhs
          logical                     :: bool
          
          bool = (REAL(lhs%f0,PR).EQ.REAL(rhs%f0,PR))
          return    
        end function hyd_eq_hyd

        function hyd_eq_dble(lhs, rhs) result(bool)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_eq_dble" :: hyd_eq_dble

          
          TYPE(hyperdual), intent(in) :: lhs
          real(PR), intent(in)        :: rhs
          logical :: bool
          
          bool = (REAL(lhs%f0,PR).EQ.rhs)
        
        end function hyd_eq_dble

        function dble_eq_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_eq_hyd" :: dble_eq_hyd

          
          real(PR), intent(in)        :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          logical                     :: bool
          
          bool = (REAL(rhs%f0,PR).EQ.lhs)
        
        end function dble_eq_hyd


      !----- NE operator (/=)
        function hyd_ne_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_ne_hyd" :: hyd_ne_hyd
        
          
          TYPE(hyperdual), intent(in) :: lhs, rhs
          logical                     :: bool
                  
          bool = (.NOT.(lhs.EQ.rhs))
                      
        end function hyd_ne_hyd

        function hyd_ne_dble(lhs, rhs) result(bool)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_ne_dble" :: hyd_ne_dble

          
          TYPE(hyperdual), intent(in) :: lhs
          real(PR), intent(in)        :: rhs
          logical                     :: bool

          bool = (REAL(lhs%f0,PR).NE.rhs)

        end function hyd_ne_dble

        function dble_ne_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_ne_hyd" :: dble_ne_hyd

          
          real(PR), intent(in)        :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          logical                     :: bool

          bool = (REAL(rhs%f0,PR).NE.lhs)

        end function dble_ne_hyd


        function hyd_gt_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "hyd_gt_hyd" :: hyd_gt_hyd
          
          
          logical :: bool
          type(hyperdual), intent(in) :: lhs
          type(hyperdual), intent(in) :: rhs
          bool = lhs%f0 > rhs%f0
          return 

        end function

        function hyd_gt_dble(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "hyd_gt_dble" :: hyd_gt_dble
          
          
          logical :: bool
          double precision, intent(in) :: rhs
          type(hyperdual), intent(in) :: lhs
          bool = lhs%f0 > rhs
          return

        end function

        function dble_gt_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "dble_gt_hyd" :: dble_gt_hyd
          
          
          logical :: bool
          type(hyperdual), intent(in) :: rhs
          double precision, intent(in) :: lhs
          bool = lhs > rhs%f0
          return   
        end function


      !----- GE operator (>=)
        function hyd_ge_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "hyd_ge_hyd" :: hyd_ge_hyd
          
          
          logical :: bool
          type(hyperdual), intent(in) :: lhs
          type(hyperdual), intent(in) :: rhs
          bool = lhs%f0 >= rhs%f0 
          return   

        end function

        function hyd_ge_dble(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "hyd_ge_dble" :: hyd_ge_dble
          
          
          logical :: bool
          double precision, intent(in) :: rhs
          type(hyperdual), intent(in) :: lhs
          bool = lhs%f0 >= rhs
          return 

        end function

        function dble_ge_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "dble_ge_hyd" :: dble_ge_hyd
          
          
          logical :: bool
          type(hyperdual), intent(in) :: rhs
          double precision, intent(in) :: lhs
          bool = lhs >= rhs%f0
          return   

        end function


      !----- LT operator (<)
        function hyd_lt_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "hyd_lt_hyd" :: hyd_lt_hyd
          
          
          logical :: bool
          type(hyperdual), intent(in) :: lhs
          type(hyperdual), intent(in) :: rhs
          bool = lhs%f0 < rhs%f0
          return

        end function

        function hyd_lt_dble(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "hyd_lt_dble" :: hyd_lt_dble
          
          
          logical :: bool
          type(hyperdual), intent(in) :: lhs
          double precision, intent(in) :: rhs
          bool = lhs%f0 < rhs
          return

        end function

        function dble_lt_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "dble_lt_hyd" :: dble_lt_hyd
          
          logical :: bool
          type(hyperdual), intent(in) :: rhs
          double precision, intent(in) :: lhs
          bool =  lhs < rhs%f0 
          return   
         
        end function


      !----- LE operator (<=)
        function hyd_le_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "hyd_le_hyd" :: hyd_le_hyd
          
          logical :: bool
          type(hyperdual), intent(in) :: lhs
          type(hyperdual), intent(in) :: rhs
          bool = lhs%f0 <= rhs%f0
          return   

        end function

        function hyd_le_dble(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "hyd_le_hyd" :: hyd_le_hyd
          
          logical :: bool
          double precision, intent(in) :: rhs
          type(hyperdual), intent(in) :: lhs
          bool = lhs%f0 <= rhs
          return   
          

        end function

        function dble_le_hyd(lhs, rhs) result(bool)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "dble_le_hyd" :: dble_le_hyd
          
          logical :: bool
          type(hyperdual), intent(in) :: rhs
          double precision, intent(in) :: lhs
          bool = lhs <= rhs%f0 
          return  
        end function


      !----------------------------------------------------------------!
      !                     ARITHMETIC OPERATORS                       !
      !----------------------------------------------------------------!

      !----- Addition operator (+)
        function hyd_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_hyd" :: hyd_plus_hyd
        
          
          TYPE(hyperdual), intent(in) :: lhs, rhs
          TYPE(hyperdual) :: res
                
          res%f0 = lhs%f0 + rhs%f0  
          res%f1 = lhs%f1 + rhs%f1
          res%f2 = lhs%f2 + rhs%f2
          res%f12 = lhs%f12 + rhs%f12
          

          return       
        end function hyd_plus_hyd
 
        function hyd_plus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_hyd_array" :: hyd_plus_hyd_array
        
          
          TYPE(hyperdual), intent(in)               :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
                
          res%f0 = lhs%f0 + rhs%f0  
          res%f1 = lhs%f1 + rhs%f1
          res%f2 = lhs%f2 + rhs%f2
          res%f12 = lhs%f12 + rhs%f12
                  
        end function hyd_plus_hyd_array
               
        function hyd_plus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_hyd_matrix" :: hyd_plus_hyd_matrix
        
          
          TYPE(hyperdual), intent(in)                 :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res
                
          res%f0 = lhs%f0 + rhs%f0  
          res%f1 = lhs%f1 + rhs%f1
          res%f2 = lhs%f2 + rhs%f2
          res%f12 = lhs%f12 + rhs%f12
                  
        end function hyd_plus_hyd_matrix
               
        function hyd_plus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_hyd_tens" :: hyd_plus_hyd_tens
        
          
          TYPE(hyperdual), intent(in)                   :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual),&
          dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res
                
          res%f0 = lhs%f0 + rhs%f0  
          res%f1 = lhs%f1 + rhs%f1
          res%f2 = lhs%f2 + rhs%f2
          res%f12 = lhs%f12 + rhs%f12
                  
        end function hyd_plus_hyd_tens
               
        function hyd_plus_dble(lhs, rhs) result(res)
        !DEC$ AT    TRIbUTES DLLEXPORT, ALIAS : "hyd_plus_dble" :: hyd_plus_dble
                
                
          TYPE(hyperdual), intent(in) :: lhs
          !real(PR), intent(in)        :: rhs
          real, intent(in)            :: rhs
          TYPE(hyperdual)             :: res
                
          res%f0 = lhs%f0 + rhs 
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12
                    
        end function hyd_plus_dble
               
        function hyd_plus_dble_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_dble_array" :: hyd_plus_dble_array
                
                
          TYPE(hyperdual), intent(in)         :: lhs
          !real(PR), dimension(:), intent(in)  :: rhs
          real, dimension(:), intent(in)  :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
                
          res%f0 = lhs%f0 + rhs 
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12
                    
        end function hyd_plus_dble_array
                
        function hyd_plus_dble_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_dble_matrix" :: hyd_plus_dble_matrix
                
                
          TYPE(hyperdual), intent(in)           :: lhs
          !real(PR), dimension(:,:), intent(in)  :: rhs
          real, dimension(:,:), intent(in)  :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res
                
          res%f0 = lhs%f0 + rhs 
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12 
                    
        end function hyd_plus_dble_matrix
                
        function hyd_plus_dble_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_dble_tens" :: hyd_plus_dble_tens
                
                
          TYPE(hyperdual), intent(in)             :: lhs
          !real(PR), dimension(:,:,:), intent(in)  :: rhs
          real, dimension(:,:,:), intent(in)  :: rhs
          TYPE(hyperdual),  dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res
          res%f0 = lhs%f0 + rhs 
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12
                    
        end function hyd_plus_dble_tens
               
        function hyd_plus_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_int" :: hyd_plus_int
                
                
          TYPE(hyperdual), intent(in) :: lhs
          integer, intent(in)         :: rhs
          TYPE(hyperdual)             :: res
                
          res%f0 = lhs%f0 + REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12     

          return
        end function hyd_plus_int
                
        function hyd_plus_int_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_int_array" :: hyd_plus_int_array
                
                
          TYPE(hyperdual), intent(in)       :: lhs
          integer, dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
                
          res%f0 = lhs%f0 + REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12      

          return
        end function hyd_plus_int_array
                
        function hyd_plus_int_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_int_matrix" :: hyd_plus_int_matrix
                
                
          TYPE(hyperdual), intent(in)         :: lhs
          integer, dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1), size(rhs,2))  :: res
                
          res%f0 = lhs%f0 + REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_plus_int_matrix
         
        function hyd_plus_int_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_plus_int_tens" :: hyd_plus_int_tens
                
                
          TYPE(hyperdual), intent(in)           :: lhs
          integer, dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
          dimension(size(rhs,1), size(rhs,2), size(rhs,3))  :: res
                
          res%f0 = lhs%f0 + REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12
          return
        end function hyd_plus_int_tens
                
        function hyd_array_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_plus_hyd" :: hyd_array_plus_hyd
                
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)               :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
                
          res%f0 = lhs%f0 + rhs%f0
          res%f1 = lhs%f1 + rhs%f1
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_array_plus_hyd
                
        function hyd_array_plus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_plus_hyd_array" :: hyd_array_plus_hyd_array
                
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs, rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
                
          res%f0 = lhs%f0 + rhs%f0
          res%f1 = lhs%f1 + rhs%f1
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return          
        end function hyd_array_plus_hyd_array
                
        function hyd_array_plus_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_plus_dble" :: hyd_array_plus_dble
                
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          real(PR), intent(in)                      :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
                
          res%f0 = lhs%f0 + rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return      
        end function hyd_array_plus_dble
                
        function hyd_array_plus_dble_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_plus_dble_array" :: hyd_array_plus_dble_array
                
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          real(PR), dimension(:), intent(in)        :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
                
          res%f0 = lhs%f0 + rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return      
        end function hyd_array_plus_dble_array
                
        function hyd_array_plus_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_plus_int" :: hyd_array_plus_int
                
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          integer, intent(in)                       :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
                
          res%f0 = lhs%f0 + REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return      
        end function hyd_array_plus_int
                
        function hyd_array_plus_int_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_plus_int_array" :: hyd_array_plus_int_array
                
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          integer, dimension(:), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
                
          res%f0 = lhs%f0 + REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return      
        end function hyd_array_plus_int_array
                
        function hyd_matrix_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_plus_hyd" :: hyd_matrix_plus_hyd

          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)                 :: rhs
          TYPE(hyperdual), dimension(size(lhs,1), size(lhs,2))  :: res

          res%f0 = rhs%f0 + lhs%f0
          res%f1 = rhs%f1 + lhs%f1 
          res%f2 = rhs%f2 + lhs%f2
          res%f12 = rhs%f12 + lhs%f12

          return
        end function hyd_matrix_plus_hyd
  
        function hyd_matrix_plus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_plus_hyd_matrix" :: hyd_matrix_plus_hyd_matrix

          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1), size(rhs,2))  :: res

          res%f0 = rhs%f0 + lhs%f0
          res%f1 = rhs%f1 + lhs%f1 
          res%f2 = rhs%f2 + lhs%f2
          res%f12 = rhs%f12 + lhs%f12


          return
        end function hyd_matrix_plus_hyd_matrix
  
        function hyd_matrix_plus_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_plus_dble" :: hyd_matrix_plus_dble

          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          real(PR), intent(in)                        :: rhs
          TYPE(hyperdual), dimension(size(lhs,1), size(lhs,2))  :: res

          res%f0 = lhs%f0 + rhs
          res%f1 = lhs%f1
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_matrix_plus_dble

        function hyd_matrix_plus_dble_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_plus_dble_matrix" :: hyd_matrix_plus_dble_matrix

          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          real(PR), dimension(:,:), intent(in)        :: rhs
          TYPE(hyperdual), dimension(size(lhs,1), size(lhs,2))  :: res

          res%f0 = lhs%f0 + rhs
          res%f1 = lhs%f1
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_matrix_plus_dble_matrix
                
        function hyd_matrix_plus_int(lhs, rhs) result(res)
        !DEC$ AT    TRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_plus_int" :: hyd_matrix_plus_int
                
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          integer, intent(in)                         :: rhs
          TYPE(hyperdual), dimension(size(lhs,1), size(lhs,2))  :: res
                
          res%f0 = lhs%f0 + REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_matrix_plus_int
                
        function hyd_matrix_plus_int_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_plus_int_matrix" :: hyd_matrix_plus_int_matrix
                
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          integer, dimension(:,:), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(lhs,1), size(lhs,2))  :: res
                
          res%f0 = lhs%f0 + REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_matrix_plus_int_matrix

        function hyd_tens_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_plus_hyd" :: hyd_tens_plus_hyd

          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)                   :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3))  :: res

          res%f0 = lhs%f0 + rhs%f0
          res%f1 = lhs%f1 + rhs%f1 
          res%f2 = lhs%f2 + rhs%f2
          res%f12 = lhs%f12 + rhs%f12

          return
        
        end function

        function hyd_tens_plus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_plus_hyd_tens" :: hyd_tens_plus_hyd_tens

          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1),size(rhs,2),size(rhs,3))  :: res

          res%f0 = lhs%f0 + rhs%f0
          res%f1 = lhs%f1 + rhs%f1 
          res%f2 = lhs%f2 + rhs%f2
          res%f12 = lhs%f12 + rhs%f12

          return
        end function hyd_tens_plus_hyd_tens

        function hyd_tens_plus_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_plus_dble" :: hyd_tens_plus_dble

          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          real(PR), intent(in)                          :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3))  :: res

          res%f0 = lhs%f0 + rhs
          res%f1 = lhs%f1
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return

        end function hyd_tens_plus_dble

        function hyd_tens_plus_dble_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_plus_dble_tens" :: hyd_tens_plus_dble_tens

          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          real(PR), dimension(:,:,:), intent(in)        :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3))  :: res

          res%f0 = lhs%f0 + rhs
          res%f1 = lhs%f1
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return

        end function hyd_tens_plus_dble_tens
                
        function hyd_tens_plus_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_plus_int" :: hyd_tens_plus_int
                
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          integer, intent(in)                           :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1), size(lhs,2), size(lhs,3))  :: res
                
          res%f0 = lhs%f0 + REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_tens_plus_int
                
        function hyd_tens_plus_int_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_plus_int_tens" :: hyd_tens_plus_int_tens
                
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          integer, dimension(:,:,:), intent(in)         :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1), size(lhs,2), size(lhs,3))  :: res
                
          res%f0 = lhs%f0 + REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_tens_plus_int_tens

        function dble_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_plus_hyd" :: dble_plus_hyd
                
          
          real(PR), intent(in)        :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: res
                
          res%f0 = rhs%f0 + lhs
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12
    
          return
        end function dble_plus_hyd
                
        function dble_plus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_plus_hyd_array" :: dble_plus_hyd_array
                
          
          real(PR), intent(in)                      :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
                
          res%f0 = rhs%f0 + lhs
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return
        end function dble_plus_hyd_array

        function dble_plus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_plus_hyd_matrix" :: dble_plus_hyd_matrix

          
          real(PR), intent(in)                        :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1), size(rhs,2))  :: res

          res%f0 = rhs%f0 + lhs
          res%f1 = rhs%f1
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return
        end function dble_plus_hyd_matrix

        function dble_plus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_plus_hyd_tens" :: dble_plus_hyd_tens
                
          
          real(PR), intent(in)                          :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res
                
          res%f0 = rhs%f0 + lhs
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return      
        end function dble_plus_hyd_tens

        function dble_array_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_array_plus_hyd" :: dble_array_plus_hyd
                
          
          real(PR), dimension(:), intent(in)        :: lhs
          TYPE(hyperdual), intent(in)               :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
                
          res%f0 = rhs%f0 + lhs
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return     
      
        end function dble_array_plus_hyd

        function dble_array_plus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_array_plus_hyd_array" :: dble_array_plus_hyd_array
                
          
          real(PR), dimension(:), intent(in)        :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
                
          res%f0 = rhs%f0 + lhs
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return
        end function dble_array_plus_hyd_array

        function dble_matrix_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_matrix_plus_hyd" :: dble_matrix_plus_hyd
                
          
          real(PR), dimension(:,:), intent(in)      :: lhs
          TYPE(hyperdual), intent(in)               :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
                
          res%f0 = rhs%f0 + lhs
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return      
        end function dble_matrix_plus_hyd

        function dble_matrix_plus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_matrix_plus_hyd_matrix" :: dble_matrix_plus_hyd_matrix
                
          
          real(PR), dimension(:,:), intent(in)        :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
                
          res%f0 = rhs%f0 + lhs
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return      
        end function dble_matrix_plus_hyd_matrix

        function dble_tens_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_tens_plus_hyd" :: dble_tens_plus_hyd
                
          
          real(PR), dimension(:,:,:), intent(in)  :: lhs
          TYPE(hyperdual), intent(in)             :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
                
          res%f0 = rhs%f0 + lhs
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return      
        end function dble_tens_plus_hyd

        function dble_tens_plus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_tens_plus_hyd_tens" :: dble_tens_plus_hyd_tens
                
          
          real(PR), dimension(:,:,:), intent(in)        :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
                
          res%f0 = rhs%f0 + lhs
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return
        end function dble_tens_plus_hyd_tens

        function int_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_plus_hyd" :: int_plus_hyd
                
          
          integer, intent(in)         :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: res
                
          res%f0 = rhs%f0 + REAL(lhs,PR)
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return    
        end function int_plus_hyd

        function int_plus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_plus_hyd_array" :: int_plus_hyd_array
                
          
          integer, intent(in)                       :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
                
          res%f0 = rhs%f0 + REAL(lhs,PR)
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12
          return    
        end function int_plus_hyd_array

        function int_plus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_plus_hyd_matrix" :: int_plus_hyd_matrix
                
          
          integer, intent(in)                         :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1), size(rhs,2))  :: res
                
          res%f0 = rhs%f0 + REAL(lhs,PR)
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return    
        end function int_plus_hyd_matrix

        function int_plus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_plus_hyd_tens" :: int_plus_hyd_tens
                
          
          integer, intent(in)                           :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1), size(rhs,2), size(rhs,3)) :: res
                
          res%f0 = rhs%f0 + REAL(lhs,PR)
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return    
        end function int_plus_hyd_tens

        function int_array_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_array_plus_hyd" :: int_array_plus_hyd
                
          
          integer, dimension(:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)       :: rhs
          TYPE(hyperdual), dimension(size(lhs)) :: res
                
          res%f0 = rhs%f0 + REAL(lhs,PR)
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return    
        end function int_array_plus_hyd

        function int_array_plus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_array_plus_hyd_array" :: int_array_plus_hyd_array
                
          
          integer, dimension(:), intent(in)         :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
                
          res%f0 = rhs%f0 + REAL(lhs,PR)
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return    
        end function int_array_plus_hyd_array

        function int_matrix_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_matrix_plus_hyd" :: int_matrix_plus_hyd
                
          
          integer, dimension(:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
                
          res%f0 = rhs%f0 + REAL(lhs,PR)
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return 
        end function int_matrix_plus_hyd

        function int_matrix_plus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_matrix_plus_hyd_matrix" :: int_matrix_plus_hyd_matrix
                
          
          integer, dimension(:,:), intent(in)         :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
                
          res%f0 = rhs%f0 + REAL(lhs,PR)
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return   
        end function int_matrix_plus_hyd_matrix

        function int_tens_plus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_tens_plus_hyd" :: int_tens_plus_hyd
                
          
          integer, dimension(:,:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)           :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
                
          res%f0 = rhs%f0 + REAL(lhs,PR)
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return    
        end function int_tens_plus_hyd

        function int_tens_plus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_tens_plus_hyd_tens" :: int_tens_plus_hyd_tens
                
          
          integer, dimension(:,:,:), intent(in)         :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
                
          res%f0 = rhs%f0 + REAL(lhs,PR)
          res%f1 = rhs%f1 
          res%f2 = rhs%f2 
          res%f12 = rhs%f12

          return    
        end function int_tens_plus_hyd_tens


      !----- Subtraction operator (-)
        function hyd_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_hyd" :: hyd_minus_hyd
        
          
          TYPE(hyperdual), intent(in) :: lhs, rhs
          TYPE(hyperdual)             :: res
          
          res%f0 = lhs%f0 - rhs%f0
          res%f1 = lhs%f1 - rhs%f1
          res%f2 = lhs%f2  - rhs%f2
          res%f12 = lhs%f12 - rhs%f12

          return
        end function hyd_minus_hyd
        
        function hyd_minus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_hyd_array" :: hyd_minus_hyd_array
        
          
          TYPE(hyperdual), intent(in)               :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
          
          res%f0 = lhs%f0 - rhs%f0
          res%f1 = lhs%f1 - rhs%f1
          res%f2 = lhs%f2  - rhs%f2
          res%f12 = lhs%f12 - rhs%f12

          return
            
        end function hyd_minus_hyd_array
        
        function hyd_minus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_hyd_matrix" :: hyd_minus_hyd_matrix
        
          
          TYPE(hyperdual), intent(in)                 :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res
          
          res%f0 = lhs%f0 - rhs%f0
          res%f1 = lhs%f1 - rhs%f1
          res%f2 = lhs%f2  - rhs%f2
          res%f12 = lhs%f12 - rhs%f12

          return
        end function hyd_minus_hyd_matrix
        
        function hyd_minus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_hyd_tens" :: hyd_minus_hyd_tens
        
          
          TYPE(hyperdual), intent(in)                   :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res
          
          res%f0 = lhs%f0 - rhs%f0
          res%f1 = lhs%f1 - rhs%f1
          res%f2 = lhs%f2  - rhs%f2
          res%f12 = lhs%f12 - rhs%f12

          return
            
        end function hyd_minus_hyd_tens
        
        function hyd_minus_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_dble" :: hyd_minus_dble
        
          
          TYPE(hyperdual), intent(in) :: lhs
          real(PR), intent(in)        :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = lhs%f0 - rhs
          res%f1 = lhs%f1
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
            
        end function hyd_minus_dble
        
        function hyd_minus_dble_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_dble_array" :: hyd_minus_dble_array
        
          
          TYPE(hyperdual), intent(in)               :: lhs
          real(PR), dimension(:), intent(in)        :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
        
          res%f0 = lhs%f0 - rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_minus_dble_array

     
        
        
        function hyd_minus_dble_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_dble_matrix" :: hyd_minus_dble_matrix
        
          
          TYPE(hyperdual), intent(in)                 :: lhs
          real(PR), dimension(:,:), intent(in)        :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res
        
          res%f0 = lhs%f0 - rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_minus_dble_matrix
        
        function hyd_minus_dble_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_dble_tens" :: hyd_minus_dble_tens
        
          
          TYPE(hyperdual), intent(in)                   :: lhs
          real(PR), dimension(:,:,:), intent(in)        :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res
        
          res%f0 = lhs%f0 - rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_minus_dble_tens
                
        function hyd_minus_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_int" :: hyd_minus_int
        
          
          TYPE(hyperdual), intent(in) :: lhs
          integer, intent(in)         :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = lhs%f0 - REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_minus_int
                 
        function hyd_minus_int_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_int_array" :: hyd_minus_int_array
        
          
          TYPE(hyperdual), intent(in)       :: lhs
          integer, dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
        
          res%f0 = lhs%f0 - REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_minus_int_array               

        function hyd_minus_int_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_int_matrix" :: hyd_minus_int_matrix
        
          
          TYPE(hyperdual), intent(in)         :: lhs
          integer, dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1), size(rhs,2)) :: res
        
          res%f0 = lhs%f0 - REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_minus_int_matrix
        
        
    
                
        function hyd_minus_int_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_minus_int_tens" :: hyd_minus_int_tens
        
          
          TYPE(hyperdual), intent(in)           :: lhs
          integer, dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1), size(rhs,2), size(rhs,3))   :: res
        
          res%f0 = lhs%f0 - REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_minus_int_tens

        function hyd_array_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_minus_hyd" :: hyd_array_minus_hyd
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)               :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
        
          res%f0 = lhs%f0 - rhs%f0
          res%f1 = lhs%f1 - rhs%f1
          res%f2 = lhs%f2 - rhs%f2
          res%f12 = lhs%f12 - rhs%f12

          return
        end function hyd_array_minus_hyd
        
        function hyd_array_minus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_minus_hyd_array" :: hyd_array_minus_hyd_array
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs, rhs
          TYPE(hyperdual), dimension(size(lhs)) :: res
        
          res%f0 = lhs%f0 - rhs%f0
          res%f1 = lhs%f1 - rhs%f1
          res%f2 = lhs%f2 - rhs%f2
          res%f12 = lhs%f12 - rhs%f12

          return
        end function hyd_array_minus_hyd_array
        
        function hyd_array_minus_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_minus_dble" :: hyd_array_minus_dble
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          real(PR), intent(in)                      :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
        
          res%f0 = lhs%f0 - rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_array_minus_dble
        
        function hyd_array_minus_dble_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_minus_dble_array" :: hyd_array_minus_dble_array
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          real(PR), dimension(:), intent(in)        :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
        
          res%f0 = lhs%f0 - rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_array_minus_dble_array
        
        function hyd_array_minus_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_minus_int" :: hyd_array_minus_int
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          integer, intent(in)                       :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
        
          res%f0 = lhs%f0 - REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_array_minus_int
        
        function hyd_array_minus_int_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_minus_int_array" :: hyd_array_minus_int_array
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          integer, dimension(:), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
        
          res%f0 = lhs%f0 - REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12

          return
        end function hyd_array_minus_int_array

        function hyd_matrix_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_minus_hyd" :: hyd_matrix_minus_hyd
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)                 :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = lhs%f0 - rhs%f0
          res%f1 = lhs%f1 - rhs%f1
          res%f2 = lhs%f2 - rhs%f2
          res%f12 = lhs%f12 - rhs%f12

          return
        end function hyd_matrix_minus_hyd
        
        function hyd_matrix_minus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_minus_hyd_matrix" :: hyd_matrix_minus_hyd_matrix
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs, rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = lhs%f0 - rhs%f0
          res%f1 = lhs%f1 - rhs%f1
          res%f2 = lhs%f2 - rhs%f2
          res%f12 = lhs%f12 - rhs%f12

          return
          
        end function hyd_matrix_minus_hyd_matrix
        
        function hyd_matrix_minus_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_minus_dble" :: hyd_matrix_minus_dble
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          real(PR), intent(in)                        :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = lhs%f0 - rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2 
          res%f12 = lhs%f12 

          return
        end function hyd_matrix_minus_dble
        
        function hyd_matrix_minus_dble_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_minus_dble_matrix" :: hyd_matrix_minus_dble_matrix
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          real(PR), dimension(:,:), intent(in)        :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = lhs%f0 - rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2
          res%f12 = lhs%f12 

          return
        end function hyd_matrix_minus_dble_matrix
        
        function hyd_matrix_minus_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_minus_int" :: hyd_matrix_minus_int
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          integer, intent(in)                         :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = lhs%f0 - REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2
          res%f12 = lhs%f12 

          return
        end function hyd_matrix_minus_int
        
        function hyd_matrix_minus_int_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_minus_int_matrix" :: hyd_matrix_minus_int_matrix
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          integer, dimension(:,:), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = lhs%f0 - REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2
          res%f12 = lhs%f12 

          return
        end function hyd_matrix_minus_int_matrix

        function hyd_tens_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_minus_hyd" :: hyd_tens_minus_hyd
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)                   :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
        
          res%f0 = lhs%f0 - rhs%f0
          res%f1 = lhs%f1 - rhs%f1
          res%f2 = lhs%f2 - rhs%f2
          res%f12 = lhs%f12 - rhs%f12

          return
        end function hyd_tens_minus_hyd
        
        function hyd_tens_minus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_minus_hyd_tens" :: hyd_tens_minus_hyd_tens
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs, rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
        
          res%f0 = lhs%f0 - rhs%f0
          res%f1 = lhs%f1 - rhs%f1
          res%f2 = lhs%f2 - rhs%f2
          res%f12 = lhs%f12 - rhs%f12

          return
        end function hyd_tens_minus_hyd_tens
        
        function hyd_tens_minus_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_minus_dble" :: hyd_tens_minus_dble
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          real(PR), intent(in)                          :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
        
          res%f0 = lhs%f0 - rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2
          res%f12 = lhs%f12 

          return
        end function hyd_tens_minus_dble
        
        function hyd_tens_minus_dble_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_minus_dble_tens" :: hyd_tens_minus_dble_tens
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          real(PR), dimension(:,:,:), intent(in)        :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
        
          res%f0 = lhs%f0 - rhs
          res%f1 = lhs%f1 
          res%f2 = lhs%f2
          res%f12 = lhs%f12 

          return
        end function hyd_tens_minus_dble_tens
        
        function hyd_tens_minus_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_minus_int" :: hyd_tens_minus_int
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          integer, intent(in)                           :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
        
          res%f0 = lhs%f0 - REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2
          res%f12 = lhs%f12 

          return
        end function hyd_tens_minus_int
        
        function hyd_tens_minus_int_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_minus_int_tens" :: hyd_tens_minus_int_tens
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          integer, dimension(:,:,:), intent(in)         :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
        
          res%f0 = lhs%f0 - REAL(rhs,PR)
          res%f1 = lhs%f1 
          res%f2 = lhs%f2
          res%f12 = lhs%f12 

          return
        end function hyd_tens_minus_int_tens

        function dble_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_minus_hyd" :: dble_minus_hyd
        
          
          real(PR), intent(in)        :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = lhs - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12

          return
        end function dble_minus_hyd
        
        function dble_minus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_minus_hyd_array" :: dble_minus_hyd_array
        
          real(PR), intent(in)                      :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res

          res%f0 = lhs - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12 

          return
        end function dble_minus_hyd_array
        
        function dble_minus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_minus_hyd_matrix" :: dble_minus_hyd_matrix
        
          real(PR), intent(in)                        :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res

          res%f0 = lhs - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12
          return
        end function dble_minus_hyd_matrix
        
        function dble_minus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_minus_hyd_tens" :: dble_minus_hyd_tens
        
          real(PR), intent(in)                          :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res

          res%f0 = lhs - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12
          return
        end function dble_minus_hyd_tens
        
        function dble_array_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_array_minus_hyd" :: dble_array_minus_hyd
        
          real(PR), dimension(:), intent(in)  :: lhs
          TYPE(hyperdual), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(lhs)) :: res

          res%f0 = lhs - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12

          return
        end function dble_array_minus_hyd
        
        function dble_array_minus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_array_minus_hyd_array" :: dble_array_minus_hyd_array
        
          real(PR), dimension(:), intent(in)        :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res

          res%f0 = lhs - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12

          return
        end function dble_array_minus_hyd_array
        
        function dble_matrix_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_matrix_minus_hyd" :: dble_matrix_minus_hyd
        
          real(PR), dimension(:,:), intent(in)  :: lhs
          TYPE(hyperdual), intent(in)           :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res

          res%f0 = lhs - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12

          return
        end function dble_matrix_minus_hyd
        
        function dble_matrix_minus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_matrix_minus_hyd_matrix" :: dble_matrix_minus_hyd_matrix
        
          real(PR), dimension(:,:), intent(in)        :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res

          res%f0 = lhs - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12

          return
        end function dble_matrix_minus_hyd_matrix
        
        function dble_tens_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_tens_minus_hyd" :: dble_tens_minus_hyd
        
          real(PR), dimension(:,:,:), intent(in)  :: lhs
          TYPE(hyperdual), intent(in)             :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res

          res%f0 = lhs - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12

          return
        end function dble_tens_minus_hyd

        function dble_tens_minus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_tens_minus_hyd_tens" :: dble_tens_minus_hyd_tens
        
          real(PR), dimension(:,:,:), intent(in)        :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res

          res%f0 = lhs - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12

          return
        end function dble_tens_minus_hyd_tens
        
        function int_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_minus_hyd" :: int_minus_hyd
        
          
          integer, intent(in)         :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = REAL(lhs,PR) - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12 

          return
        end function int_minus_hyd
        
        function int_minus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_minus_hyd_array" :: int_minus_hyd_array
        
          
          integer, intent(in)                       :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
        
          res%f0 = REAL(lhs,PR) - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12 

          return
        end function int_minus_hyd_array

        function int_minus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_minus_hyd_matrix" :: int_minus_hyd_matrix
        
          
          integer, intent(in)                         :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1), size(rhs,2))  :: res
        
          res%f0 = REAL(lhs,PR) - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12 

          return
        end function int_minus_hyd_matrix

        function int_minus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_minus_hyd_tens" :: int_minus_hyd_tens
        
          
          integer, intent(in)                           :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1), size(rhs,2), size(rhs,3))  :: res
        
          res%f0 = REAL(lhs,PR) - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12
          return
        end function int_minus_hyd_tens
        
        function int_array_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_array_minus_hyd" :: int_array_minus_hyd
        
          
          integer, dimension(:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)       :: rhs
          TYPE(hyperdual), dimension(size(lhs)) :: res
        
          res%f0 = REAL(lhs,PR) - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12

          return
        end function int_array_minus_hyd
        
        function int_array_minus_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_array_minus_hyd_array" :: int_array_minus_hyd_array
        
          
          integer, dimension(:), intent(in)         :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
        
          res%f0 = REAL(lhs,PR) - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12

          return
        end function int_array_minus_hyd_array
        
        function int_matrix_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_matrix_minus_hyd" :: int_matrix_minus_hyd
        
          
          integer, dimension(:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(lhs,1), size(lhs,2))  :: res
        
          res%f0 = REAL(lhs,PR) - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12
          return
        end function int_matrix_minus_hyd
        
        function int_matrix_minus_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_matrix_minus_hyd_matrix" :: int_matrix_minus_hyd_matrix
        
          
          integer, dimension(:,:), intent(in)         :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(lhs,1), size(lhs,2))  :: res
        
          res%f0 = REAL(lhs,PR) - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12
          return
        end function int_matrix_minus_hyd_matrix
        
        function int_tens_minus_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_tens_minus_hyd" :: int_tens_minus_hyd
        
          
          integer, dimension(:,:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)           :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1), size(lhs,2), size(lhs,3))  :: res
        
          res%f0 = REAL(lhs,PR) - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12

          return
        end function int_tens_minus_hyd
        
        function int_tens_minus_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_tens_minus_hyd_tens" :: int_tens_minus_hyd_tens
        
          
          integer, dimension(:,:,:), intent(in)         :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1), size(lhs,2), size(lhs,3))  :: res
        
          res%f0 = REAL(lhs,PR) - rhs%f0 
          res%f1 = rhs%f1 
          res%f2 = rhs%f2
          res%f12 = rhs%f12 

          return
        end function int_tens_minus_hyd_tens

        function minus_hyd(rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "minus_hyd" :: minus_hyd
        
          
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = - rhs%f0 
          res%f1 = - rhs%f1
          res%f2 = - rhs%f2
          res%f12 = - rhs%f12 

          return
        end function minus_hyd
        
        function minus_hyd_array(rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "minus_hyd_array" :: minus_hyd_array
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
        
          res%f0 = - rhs%f0 
          res%f1 = - rhs%f1
          res%f2 = - rhs%f2
          res%f12 = - rhs%f12 

          
            
          return
        end function minus_hyd_array
        
        function minus_hyd_matrix(rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "minus_hyd_matrix" :: minus_hyd_matrix
            
 
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1), size(rhs,2))  :: res
        
          res%f0 = - rhs%f0 
          res%f1 = - rhs%f1
          res%f2 = - rhs%f2
          res%f12 = -  rhs%f12 

          return
        end function minus_hyd_matrix
        
        function minus_hyd_tens(rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "minus_hyd_tens" :: minus_hyd_tens
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1),size(rhs,2),size(rhs,3))  :: res
        
          res%f0 = - rhs%f0 
          res%f1 = - rhs%f1
          res%f2 = - rhs%f2
          res%f12 = -  rhs%f12 

          return
        end function minus_hyd_tens


      !----- Multiplication operator (*)
        function hyd_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_hyd" :: hyd_mul_hyd
        
          
          TYPE(hyperdual), intent(in) :: lhs, rhs
          TYPE(hyperdual)             :: res
        
          !res%f0 = lhs%f0 * rhs%f0 - lhs%f1 * rhs%f1
          !res%f1 = lhs%f0 * rhs%f1 + lhs%f1 * rhs%f0
          res%f0 = lhs%f0 * rhs%f0
          res%f1 = (lhs%f0 * rhs%f1) + (lhs%f1 * rhs%f0)
          res%f2 = (lhs%f0 * rhs%f2) + (lhs%f2 * rhs%f0)
          res%f12 = (lhs%f0 * rhs%f12) + (lhs%f1 * rhs%f2) + (lhs%f2 * rhs%f1) + (lhs%f12 * rhs%f0)
         
         return 
        end function hyd_mul_hyd
        
        function hyd_mul_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_hyd_array" :: hyd_mul_hyd_array
        
          
          TYPE(hyperdual), intent(in)               :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
        
          res%f0 = lhs%f0 * rhs%f0
          res%f1 = (lhs%f0 * rhs%f1) + (lhs%f1 * rhs%f0)
          res%f2 = (lhs%f0 * rhs%f2) + (lhs%f2 * rhs%f0)
          res%f12 = (lhs%f0 * rhs%f12) + (lhs%f1 * rhs%f2) + (lhs%f2 * rhs%f1) + (lhs%f12 * rhs%f0)
          return
        end function hyd_mul_hyd_array
        
        function hyd_mul_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_hyd_matrix" :: hyd_mul_hyd_matrix
        
          
          TYPE(hyperdual), intent(in)                 :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res
        
          res%f0 = lhs%f0 * rhs%f0
          res%f1 = (lhs%f0 * rhs%f1) + (lhs%f1 * rhs%f0)
          res%f2 = (lhs%f0 * rhs%f2) + (lhs%f2 * rhs%f0)
          res%f12 = (lhs%f0 * rhs%f12) + (lhs%f1 * rhs%f2) + (lhs%f2 * rhs%f1) + (lhs%f12 * rhs%f0)
          return
        end function hyd_mul_hyd_matrix
        
        function hyd_mul_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_hyd_tens" :: hyd_mul_hyd_tens
        
          
          TYPE(hyperdual), intent(in)                   :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), &
          dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res
          res%f0 = lhs%f0 * rhs%f0
          res%f1 = (lhs%f0 * rhs%f1) + (lhs%f1 * rhs%f0)
          res%f2 = (lhs%f0 * rhs%f2) + (lhs%f2 * rhs%f0)
          res%f12 = (lhs%f0 * rhs%f12) + (lhs%f1 * rhs%f2) + (lhs%f2 * rhs%f1) + (lhs%f12 * rhs%f0)
          return
        end function hyd_mul_hyd_tens

        function hyd_mul_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_dble" :: hyd_mul_dble
        
          
          TYPE(hyperdual), intent(in) :: lhs
          real(PR), intent(in)        :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = lhs%f0 * rhs
          res%f1 = lhs%f1 * rhs
          res%f2 = lhs%f2 * rhs
          res%f12 = lhs%f12 * rhs
          return
        end function hyd_mul_dble
        
        function hyd_mul_dble_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_dble_array" :: hyd_mul_dble_array
        
          
          TYPE(hyperdual), intent(in)         :: lhs
          real(PR), dimension(:), intent(in)  :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
        
          res%f0 = lhs%f0 * rhs
          res%f1 = lhs%f1 * rhs
          res%f2 = lhs%f2 * rhs
          res%f12 = lhs%f12 * rhs
            
          return
        end function hyd_mul_dble_array

        function hyd_mul_dble_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_dble_matrix" :: hyd_mul_dble_matrix
        
          
          TYPE(hyperdual), intent(in)           :: lhs
          real(PR), dimension(:,:), intent(in)  :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res
        
          res%f0 = lhs%f0 * rhs
          res%f1 = lhs%f1 * rhs
          res%f2 = lhs%f2 * rhs
          res%f12 = lhs%f12 * rhs
        
          return
        end function hyd_mul_dble_matrix

        function hyd_mul_dble_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_dble_tens" :: hyd_mul_dble_tens
        
          
          TYPE(hyperdual), intent(in)             :: lhs
          real(PR), dimension(:,:,:), intent(in)  :: rhs
          TYPE(hyperdual), &
          dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res
        
          res%f0 = lhs%f0 * rhs
          res%f1 = lhs%f1 * rhs
          res%f2 = lhs%f2 * rhs
          res%f12 = lhs%f12 * rhs
        
          return
        end function hyd_mul_dble_tens

        function hyd_mul_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_int" :: hyd_mul_int
        
          
          TYPE(hyperdual), intent(in) :: lhs
          integer, intent(in)         :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = lhs%f0 * REAL(rhs,PR)
          res%f1 = lhs%f1 * REAL(rhs,PR)
          res%f2 = lhs%f2 * REAL(rhs,PR)
          res%f12 = lhs%f12 * REAL(rhs,PR)
        
          return
        end function hyd_mul_int

        function hyd_mul_int_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_int_array" :: hyd_mul_int_array
        
          
          TYPE(hyperdual), intent(in)               :: lhs
          integer, dimension(:), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(rhs))     :: res
        
          res%f0 = lhs%f0 * REAL(rhs,PR)
          res%f1 = lhs%f1 * REAL(rhs,PR)
          res%f2 = lhs%f2 * REAL(rhs,PR)
          res%f12 = lhs%f12 * REAL(rhs,PR)
        
          return
        end function hyd_mul_int_array

        function hyd_mul_int_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_int_matrix" :: hyd_mul_int_matrix
        
          
          TYPE(hyperdual), intent(in)                 :: lhs
          integer, dimension(:,:), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res
        
          res%f0 = lhs%f0 * REAL(rhs,PR)
          res%f1 = lhs%f1 * REAL(rhs,PR)
          res%f2 = lhs%f2 * REAL(rhs,PR)
          res%f12 = lhs%f12 * REAL(rhs,PR)
        
          return
        end function hyd_mul_int_matrix

        function hyd_mul_int_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_mul_int_tens" :: hyd_mul_int_tens
        
          
          TYPE(hyperdual), intent(in)                   :: lhs
          integer, dimension(:,:,:), intent(in)         :: rhs
          TYPE(hyperdual), &
            dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res
        
          res%f0 = lhs%f0 * REAL(rhs,PR)
          res%f1 = lhs%f1 * REAL(rhs,PR)
          res%f2 = lhs%f2 * REAL(rhs,PR)
          res%f12 = lhs%f12 * REAL(rhs,PR)
        
          return
        end function hyd_mul_int_tens

        function hyd_array_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_mul_hyd" :: hyd_array_mul_hyd

          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)               :: rhs
          TYPE(hyperdual), dimension(size(lhs))     :: res
          res%f0 = lhs%f0 * rhs%f0
          res%f1 = (lhs%f0 * rhs%f1) + (lhs%f1 * rhs%f0)
          res%f2 = (lhs%f0 * rhs%f2) + (lhs%f2 * rhs%f0)
          res%f12 = (lhs%f0 * rhs%f12) + (lhs%f1 * rhs%f2) + (lhs%f2 * rhs%f1) + (lhs%f12 * rhs%f0)
          return
        end function hyd_array_mul_hyd

        function hyd_array_mul_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_mul_dble" :: hyd_array_mul_dble

          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          real(PR), intent(in)                      :: rhs
          TYPE(hyperdual), dimension(size(lhs)) :: res

          res%f0 = lhs%f0 * rhs
          res%f1 = lhs%f1 * rhs
          res%f2 = lhs%f2 * rhs
          res%f12 = lhs%f12 * rhs

          return
        end function hyd_array_mul_dble

        function hyd_array_mul_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_mul_int" :: hyd_array_mul_int
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          integer, intent(in)                       :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
        
          res%f0 = lhs%f0 * REAL(rhs,PR)
          res%f1 = lhs%f1 * REAL(rhs,PR)
          res%f2 = lhs%f2 * REAL(rhs,PR)
          res%f12 = lhs%f12 * REAL(rhs,PR)
        
          return
        end function hyd_array_mul_int

        function hyd_matrix_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_mul_hyd" :: hyd_matrix_mul_hyd
          
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)                 :: rhs
          TYPE(hyperdual), dimension(size(lhs,1), size(lhs,2))  :: res
    
          res%f0 = lhs%f0 * rhs%f0
          res%f1 = (lhs%f0 * rhs%f1) + (lhs%f1 * rhs%f0)
          res%f2 = (lhs%f0 * rhs%f2) + (lhs%f2 * rhs%f0)
          res%f12 = (lhs%f0 * rhs%f12) + (lhs%f1 * rhs%f2) + (lhs%f2 * rhs%f1) + (lhs%f12 * rhs%f0)

          return
        end function hyd_matrix_mul_hyd

        function hyd_matrix_mul_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_mul_dble" :: hyd_matrix_mul_dble

          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          real(PR), intent(in)                        :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res

          res%f0 = lhs%f0 * rhs
          res%f1 = lhs%f1 * rhs
          res%f2 = lhs%f2 * rhs
          res%f12 = lhs%f12 * rhs

          return
          
        end function hyd_matrix_mul_dble

        function hyd_matrix_mul_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_mul_int" :: hyd_matrix_mul_int
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          integer, intent(in)                         :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = lhs%f0 * REAL(rhs,PR)
          res%f1 = lhs%f1 * REAL(rhs,PR)
          res%f2 = lhs%f2 * REAL(rhs,PR)
          res%f12 = lhs%f12 * REAL(rhs,PR)

          return
        
        end function hyd_matrix_mul_int

        function hyd_tens_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_mul_hyd" :: hyd_tens_mul_hyd

          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)                   :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res

          res%f0 = lhs%f0 * rhs%f0
          res%f1 = (lhs%f0 * rhs%f1) + (lhs%f1 * rhs%f0)
          res%f2 = (lhs%f0 * rhs%f2) + (lhs%f2 * rhs%f0)
          res%f12 = (lhs%f0 * rhs%f12) + (lhs%f1 * rhs%f2) + (lhs%f2 * rhs%f1) + (lhs%f12 * rhs%f0)

          return
          
        end function

        function hyd_tens_mul_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_mul_dble" :: hyd_tens_mul_dble

          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          real(PR), intent(in)                          :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res

          res%f0 = lhs%f0 * rhs
          res%f1 = lhs%f1 * rhs
          res%f2 = lhs%f2 * rhs
          res%f12 = lhs%f12 * rhs

          return
          
        end function

        function hyd_tens_mul_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_mul_int" :: hyd_tens_mul_int
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          integer, intent(in)                           :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
        
          res%f0 = lhs%f0 * REAL(rhs,PR)
          res%f1 = lhs%f1 * REAL(rhs,PR)
          res%f2 = lhs%f2 * REAL(rhs,PR)
          res%f12 = lhs%f12 * REAL(rhs,PR)

          return
        end function hyd_tens_mul_int

        function dble_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_mul_hyd" :: dble_mul_hyd
        
          
          real(PR), intent(in)        :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = rhs%f0 * lhs
          res%f1 = rhs%f1 * lhs
          res%f2 = rhs%f2 * lhs
          res%f12 = rhs%f12 * lhs

          return
        end function dble_mul_hyd
        
        function dble_mul_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_mul_hyd_array" :: dble_mul_hyd_array
        
          
          real(PR), intent(in)                      :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
        
          res%f0 = rhs%f0 * lhs
          res%f1 = rhs%f1 * lhs
          res%f2 = rhs%f2 * lhs
          res%f12 = rhs%f12 * lhs

          return
        
        end function dble_mul_hyd_array
        
        function dble_mul_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_mul_hyd_matrix" :: dble_mul_hyd_matrix
        
          
          real(PR), intent(in)                        :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res
        
          res%f0 = rhs%f0 * lhs
          res%f1 = rhs%f1 * lhs
          res%f2 = rhs%f2 * lhs
          res%f12 = rhs%f12 * lhs

          return
        end function dble_mul_hyd_matrix
        
        function dble_mul_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_mul_hyd_tens" :: dble_mul_hyd_tens
        
          
          real(PR), intent(in)                          :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res
        
          res%f0 = rhs%f0 * lhs
          res%f1 = rhs%f1 * lhs
          res%f2 = rhs%f2 * lhs
          res%f12 = rhs%f12 * lhs

          return
        end function dble_mul_hyd_tens
        
        function dble_array_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_array_mul_hyd" :: dble_array_mul_hyd
        
          
          real(PR), dimension(:), intent(in)  :: lhs
          TYPE(hyperdual), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(lhs)) :: res
        
          res%f0 = rhs%f0 * lhs
          res%f1 = rhs%f1 * lhs
          res%f2 = rhs%f2 * lhs
          res%f12 = rhs%f12 * lhs

          return
        end function dble_array_mul_hyd
        
        function dble_matrix_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_matrix_mul_hyd" :: dble_matrix_mul_hyd
        
          
          real(PR), dimension(:,:), intent(in)  :: lhs
          TYPE(hyperdual), intent(in)           :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = rhs%f0 * lhs
          res%f1 = rhs%f1 * lhs
          res%f2 = rhs%f2 * lhs
          res%f12 = rhs%f12 * lhs

          return
        end function dble_matrix_mul_hyd
        
        function dble_tens_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_tens_mul_hyd" :: dble_tens_mul_hyd
        
          
          real(PR), dimension(:,:,:), intent(in)  :: lhs
          TYPE(hyperdual), intent(in)             :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
        
          res%f0 = rhs%f0 * lhs
          res%f1 = rhs%f1 * lhs
          res%f2 = rhs%f2 * lhs
          res%f12 = rhs%f12 * lhs

          return
        end function dble_tens_mul_hyd
        
        function int_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_mul_hyd" :: int_mul_hyd
        
          
          integer, intent(in)         :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = rhs%f0 * REAL(lhs,PR)
          res%f1 = rhs%f1 * REAL(lhs,PR)
          res%f2 = rhs%f2 * REAL(lhs,PR)
          res%f12 = rhs%f12 * REAL(lhs,PR)

          return
        end function int_mul_hyd

        function int_mul_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_mul_hyd_array" :: int_mul_hyd_array
        
          
          integer, intent(in)                       :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res
        
          res%f0 = rhs%f0 * REAL(lhs,PR)
          res%f1 = rhs%f1 * REAL(lhs,PR)
          res%f2 = rhs%f2 * REAL(lhs,PR)
          res%f12 = rhs%f12 * REAL(lhs,PR)

          return
        end function int_mul_hyd_array

        function int_mul_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_mul_hyd_matrix" :: int_mul_hyd_matrix
        
          
          integer, intent(in)                         :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2)) :: res
        
          res%f0 = rhs%f0 * REAL(lhs,PR)
          res%f1 = rhs%f1 * REAL(lhs,PR)
          res%f2 = rhs%f2 * REAL(lhs,PR)
          res%f12 = rhs%f12 * REAL(lhs,PR)
        
          return
        end function int_mul_hyd_matrix

        function int_mul_hyd_tens(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_mul_hyd_tens" :: int_mul_hyd_tens
        
          
          integer, intent(in)                           :: lhs
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: res
        
          res%f0 = rhs%f0 * REAL(lhs,PR)
          res%f1 = rhs%f1 * REAL(lhs,PR)
          res%f2 = rhs%f2 * REAL(lhs,PR)
          res%f12 = rhs%f12 * REAL(lhs,PR)

          return
        end function int_mul_hyd_tens

        function int_array_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_array_mul_hyd" :: int_array_mul_hyd
        
          
          integer, dimension(:), intent(in)       :: lhs
          TYPE(hyperdual), intent(in)             :: rhs
          TYPE(hyperdual), dimension(size(lhs)) :: res
        
          res%f0 = rhs%f0 * REAL(lhs,PR)
          res%f1 = rhs%f1 * REAL(lhs,PR)
          res%f2 = rhs%f2 * REAL(lhs,PR)
          res%f12 = rhs%f12 * REAL(lhs,PR)

          return
        end function int_array_mul_hyd

        function int_matrix_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_matrix_mul_hyd" :: int_matrix_mul_hyd
        
          
          integer, dimension(:,:), intent(in)       :: lhs
          TYPE(hyperdual), intent(in)               :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = rhs%f0 * REAL(lhs,PR)
          res%f1 = rhs%f1 * REAL(lhs,PR)
          res%f2 = rhs%f2 * REAL(lhs,PR)
          res%f12 = rhs%f12 * REAL(lhs,PR)

          return
        end function int_matrix_mul_hyd

        function int_tens_mul_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_tens_mul_hyd" :: int_tens_mul_hyd
        
          
          integer, dimension(:,:,:), intent(in)       :: lhs
          TYPE(hyperdual), intent(in)                 :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
        
          res%f0 = rhs%f0 * REAL(lhs,PR)
          res%f1 = rhs%f1 * REAL(lhs,PR)
          res%f2 = rhs%f2 * REAL(lhs,PR)
          res%f12 = rhs%f12 * REAL(lhs,PR)

          return
        end function int_tens_mul_hyd

        
      !----- Division operator (/)
        function hyd_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_div_hyd" :: hyd_div_hyd
        
          
          TYPE(hyperdual), intent(in) :: lhs, rhs
          TYPE(hyperdual)             :: res, inv
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = hyd_mul_hyd(lhs, inv)

          return
          
        end function hyd_div_hyd
        
        function hyd_div_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_div_dble" :: hyd_div_dble
        
          
          TYPE(hyperdual), intent(in) :: lhs
          real(PR), intent(in)        :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = lhs%f0 / rhs
          res%f1 = lhs%f1 / rhs
          res%f2 = lhs%f2 / rhs
          res%f12 = lhs%f12 / rhs

          return
          
        end function hyd_div_dble
        
        function hyd_div_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_div_int" :: hyd_div_int
        
           
          TYPE(hyperdual), intent(in) :: lhs
          integer, intent(in)         :: rhs
          TYPE(hyperdual)             :: res
        
          res%f0 = lhs%f0 / REAL(rhs, PR)
          res%f1 = lhs%f1 / REAL(rhs, PR)
          res%f2 = lhs%f2 / REAL(rhs, PR)
          res%f12 = lhs%f12 / REAL(rhs, PR)
          
          return  
           
        end function hyd_div_int
        
        function hyd_array_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_div_hyd" :: hyd_array_div_hyd
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)               :: rhs
          TYPE(hyperdual), dimension(size(lhs))     :: res
          TYPE(hyperdual)                           :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = hyd_array_mul_hyd(lhs, inv)
          
          return
          
        end function hyd_array_div_hyd
        
        function hyd_array_div_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_div_dble" :: hyd_array_div_dble
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          real(PR), intent(in)                      :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
        
          res%f0 = lhs%f0 / rhs
          res%f1 = lhs%f1 / rhs
          res%f2 = lhs%f2 / rhs
          res%f12 = lhs%f12 / rhs
          
          return
          
        end function hyd_array_div_dble
        
        function hyd_array_div_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_array_div_int" :: hyd_array_div_int
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          integer, intent(in)                       :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
        
          res%f0 = lhs%f0 / REAL(rhs, PR)
          res%f1 = lhs%f1 / REAL(rhs, PR)
          res%f2 = lhs%f2 / REAL(rhs, PR)
          res%f12 = lhs%f12 / REAL(rhs, PR)

          return
            
        end function hyd_array_div_int
        
        function hyd_matrix_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_div_hyd" :: hyd_matrix_div_hyd
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)                 :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
          TYPE(hyperdual)                          :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = hyd_matrix_mul_hyd(lhs, inv)
            
          return
          
        end function hyd_matrix_div_hyd
        
        function hyd_matrix_div_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_div_dble" :: hyd_matrix_div_dble
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          real(PR), intent(in)                        :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = lhs%f0 / rhs
          res%f1 = lhs%f1 / rhs
          res%f2 = lhs%f2 / rhs
          res%f12 = lhs%f12 / rhs
          
          return
          
        end function hyd_matrix_div_dble
        
        function hyd_matrix_div_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_matrix_div_int" :: hyd_matrix_div_int
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          integer, intent(in)                         :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
        
          res%f0 = lhs%f0 / REAL(rhs, PR)
          res%f1 = lhs%f1 / REAL(rhs, PR)
          res%f2 = lhs%f2 / REAL(rhs, PR)
          res%f12 = lhs%f12 / REAL(rhs, PR)

          return
          
        end function hyd_matrix_div_int
        
        function hyd_tens_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_div_hyd" :: hyd_tens_div_hyd
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)                   :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3))  :: res
          TYPE(hyperdual)                       :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = hyd_tens_mul_hyd(lhs, inv)
            
          return
          
        end function hyd_tens_div_hyd
        
        function hyd_tens_div_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_div_dble" :: hyd_tens_div_dble
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          real(PR), intent(in)                          :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3))  :: res
        
          res%f0 = lhs%f0 / rhs
          res%f1 = lhs%f1 / rhs
          res%f2 = lhs%f2 / rhs
          res%f12 = lhs%f12 / rhs
          
          return

        end function hyd_tens_div_dble
        
        function hyd_tens_div_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_tens_div_int" :: hyd_tens_div_int
        
          
          TYPE(hyperdual), dimension(:,:,:), intent(in) :: lhs
          integer, intent(in)                           :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
        
          res%f0 = lhs%f0 / REAL(rhs, PR)
          res%f1 = lhs%f1 / REAL(rhs, PR)
          res%f2 = lhs%f2 / REAL(rhs, PR)
          res%f12 = lhs%f12 / REAL(rhs, PR)

          return
          
        end function hyd_tens_div_int
        
        function dble_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_div_hyd" :: dble_div_hyd

          
          real(PR), intent(in)        :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: res
          TYPE(hyperdual)             :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = dble_mul_hyd(lhs, inv)

          return
          
        end function dble_div_hyd
        
        function int_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_div_hyd" :: int_div_hyd
        
          
          integer, intent(in)         :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: res
          TYPE(hyperdual)             :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = int_mul_hyd(lhs, inv)

          return
          
        end function int_div_hyd    

        function dble_array_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_array_div_hyd" :: dble_array_div_hyd
        
          
          real(PR), dimension(:), intent(in)      :: lhs
          TYPE(hyperdual), intent(in)             :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res
          TYPE(hyperdual)                         :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = dble_array_mul_hyd(lhs, inv)

          return
          
        end function dble_array_div_hyd
        
        function int_array_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_array_div_hyd" :: int_array_div_hyd
        
          
          integer, dimension(:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)       :: rhs
          TYPE(hyperdual), dimension(size(lhs,1)) :: res
          TYPE(hyperdual)                         :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = int_array_mul_hyd(lhs, inv)

          return
          
        end function int_array_div_hyd    

        function dble_matrix_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_matrix_div_hyd" :: dble_matrix_div_hyd
        
          
          real(PR), dimension(:,:), intent(in)  :: lhs
          TYPE(hyperdual), intent(in)           :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
          TYPE(hyperdual)                       :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = dble_matrix_mul_hyd(lhs, inv)

          return
          
        end function dble_matrix_div_hyd
        
        function int_matrix_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_matrix_div_hyd" :: int_matrix_div_hyd
        
          
          integer, dimension(:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)         :: rhs
          TYPE(hyperdual), dimension(size(lhs,1),size(lhs,2)) :: res
          TYPE(hyperdual)                     :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = int_matrix_mul_hyd(lhs, inv)
          
        
          return
        end function int_matrix_div_hyd       

        function dble_tens_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "dble_tens_div_hyd" :: dble_tens_div_hyd
        
          
          real(PR), dimension(:,:,:), intent(in)  :: lhs
          TYPE(hyperdual), intent(in)             :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
          TYPE(hyperdual)                         :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = dble_tens_mul_hyd(lhs, inv)
        
          return
          
        end function dble_tens_div_hyd
         
        function int_tens_div_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "int_tens_div_hyd" :: int_tens_div_hyd
        
          
          integer, dimension(:,:,:), intent(in) :: lhs
          TYPE(hyperdual), intent(in)           :: rhs
          TYPE(hyperdual), &
            dimension(size(lhs,1),size(lhs,2),size(lhs,3)) :: res
          TYPE(hyperdual)                       :: inv
          
          inv = hyd_pow_dble(rhs, -1.0_PR)
          res = int_tens_mul_hyd(lhs, inv)

          return
          
        end function int_tens_div_hyd       

        

      !----- POW operator (**)
        function hyd_pow_hyd(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_pow_hyd" :: hyd_pow_hyd
          
          
          TYPE(hyperdual), intent(in) :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: res

          res = exp_hyd(REAL(rhs%f0,PR)*log_hyd(lhs))

          return
        end function

        function hyd_pow_dble(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_pow_dble
          type(hyperdual) :: res
          type(hyperdual), intent(in) :: lhs
          double precision, intent(in) :: rhs
          double precision :: deriv
          double precision :: lhsval
          double precision :: tol
          lhsval = lhs%f0
          tol = exp(1.0)- 15
          if (abs(lhsval) < tol) then
               if(lhsval >= 0) then
                   lhsval = tol
               end if
               if(lhsval < 0) then
                   lhsval = -tol
               end if
          end if
          deriv = rhs * (lhsval**(rhs-1))
          res%f0 = (lhs%f0)**rhs
          res%f1 = lhs%f1 * deriv
          res%f2 = lhs%f2 * deriv
          res%f12 = lhs%f12 * deriv + rhs * (rhs-1) * lhs%f1 * lhs%f2 * (lhsval**(rhs-2))
          return
        
        end function hyd_pow_dble
        
    
        function hyd_pow_int(lhs, rhs) result(res)
        !DEC$ ATTRIbUTES DLLEXPORT, ALIAS : "hyd_pow_int" :: hyd_pow_int
          type(hyperdual) :: res
          type(hyperdual), intent(in) :: lhs
          integer, intent(in) :: rhs
          double precision :: deriv
          double precision :: lhsval
          double precision :: tol
          lhsval = lhs%f0
          tol = exp(1.0) - 15
          if (abs(lhsval) < tol) then
               if(lhsval >= 0) then
                   lhsval = tol
               end if
               if(lhsval < 0) then
                   lhsval = -tol
               end if
          end if
          deriv = rhs * (lhsval**(rhs-1))
          res%f0 = (lhs%f0)**rhs
          res%f1 = lhs%f1 * deriv
          res%f2 = lhs%f2 * deriv
          res%f12 = lhs%f12 * deriv + rhs * (rhs-1) * lhs%f1 * lhs%f2 * (lhsval**(rhs-2))
          return  
        end function hyd_pow_int
            

      !----------------------------------------------------------------!
      !                     INTRINSIC FUNCTIONS                        !
      !----------------------------------------------------------------!

      !----- DbLE (conversion to double)
        function dble_hyd(hyd_in) result(dble_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "dble_hyd" :: dble_hyd
          
          
          TYPE(hyperdual), intent(in)   :: hyd_in
          real(PR)                      :: dble_out

          dble_out = DBLE(hyd_in%f0) 
        
          return
        end function dble_hyd

        function dble_hyd_array(hyd_in) result(array_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "dble_hyd_array" :: dble_hyd_array
          
          
          TYPE(hyperdual), dimension(:), intent(in) :: hyd_in
          real(PR), dimension(size(hyd_in))           :: array_out
          integer :: i

          do i = 1, size(hyd_in)
            array_out(i) = DBLE(hyd_in(i)%f0)
          enddo
          
          return
        end function dble_hyd_array

        function dble_hyd_matrix(hyd_in) result(matrix_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "dble_hyd_matrix" :: dble_hyd_matrix
          
          
          TYPE(hyperdual), dimension(:,:), intent(in)     :: hyd_in
          real(PR), dimension(size(hyd_in,1), size(hyd_in,2)) :: X_out
          integer :: i, j

          do i = 1, size(hyd_in,1)
            do j = 1, size(hyd_in,2)
              matrix_out(i,j) = DBLE(hyd_in(i,j)%f0)
            enddo
          enddo
            
          return
        end function dble_hyd_matrix


      !----- ABS (absolute value)
        function abs_hyd(hyd_in) result(abs_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "abs_hyd" :: abs_hyd
          
          type(hyperdual), intent(in) :: hyd_in
          type(hyperdual) :: abs_out
          double precision :: n = -1.0
          double precision :: n2 = 0.0
          if (hyd_le_dble(hyd_in, n2)) then
              abs_out = hyd_mul_dble(hyd_in, n)
          else 
              abs_out = hyd_in
          end if
          return
        end function abs_hyd

        function abs_hyd_array(hyd_in) result(X_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "abs_hyd_array" :: abs_hyd_array
          
          
          TYPE(hyperdual), dimension(:), intent(in) :: hyd_in
          TYPE(hyperdual), dimension(size(hyd_in))    :: X_out
          integer :: i

          do i = 1, size(hyd_in)
            X_out(i) = abs(hyd_in(i))
          enddo
          return
        end function abs_hyd_array

        function abs_hyd_matrix(hyd_in) result(X_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "abs_hyd_matrix" :: abs_hyd_matrix
          
          
          TYPE(hyperdual), dimension(:,:), intent(in)           :: hyd_in
          TYPE(hyperdual), dimension(size(hyd_in,1), size(hyd_in,2)):: X_out
          integer :: i, j

          do i = 1, size(hyd_in,1)
            do j = 1, size(hyd_in,2)
              X_out(i,j) = abs(hyd_in(i,j))
            enddo
          enddo
          return
        end function abs_hyd_matrix


      !----- SIGN
        function sign_hyd_hyd(lhs, rhs) result(hyd_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "sign_hyd_hyd" :: sign_hyd_hyd
          
          
          TYPE(hyperdual), intent(in)   :: lhs, rhs
          TYPE(hyperdual)               :: hyd_out

          if (REAL(rhs%f0,PR).GE.(0.0_PR)) then
            hyd_out = abs(lhs)
          else
            hyd_out = -abs(lhs)
          endif
          return
        end function


      !----- MAX
        function max_hyd_hyd(lhs, rhs) result(hyd_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "max_hyd_hyd" :: max_hyd_hyd

          
          TYPE(hyperdual), intent(in) :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: hyd_out

          if (lhs.GE.rhs) then
            hyd_out = lhs
          else
            hyd_out = rhs
          endif
          
          return
          
        end function

        function max_hyd_dble(lhs, rhs) result(hyd_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "max_hyd_dble" :: max_hyd_dble

          
          TYPE(hyperdual), intent(in) :: lhs
          real(PR), intent(in)        :: rhs
          TYPE(hyperdual)             :: hyd_out

          if (lhs.GE.rhs) then
            hyd_out = lhs
          else
            hyd_out = rhs
          endif
          
          return
        end function

        function max_dble_hyd(lhs, rhs) result(hyd_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "max_dble_hyd" :: max_dble_hyd

          
          real(PR), intent(in)        :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: hyd_out

          if (lhs.GE.rhs) then
            hyd_out = lhs
          else
            hyd_out = rhs
          end if

          return
        end function


      !----- MIN
        function min_hyd_hyd(lhs, rhs) result(hyd_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "min_hyd_hyd" :: min_hyd_hyd

          
          TYPE(hyperdual), intent(in) :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: hyd_out

          if (lhs.LE.rhs) then
            hyd_out = lhs
          else
            hyd_out = rhs
          endif

          return
        end function

        function min_hyd_dble(lhs, rhs) result(hyd_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "min_hyd_dble" :: min_hyd_dble

          
          TYPE(hyperdual), intent(in) :: lhs
          real(PR), intent(in)        :: rhs
          TYPE(hyperdual)             :: hyd_out

          if (lhs.LE.rhs) then
            hyd_out = lhs
          else
            hyd_out = rhs
          endif

          return
        end function

        function min_dble_hyd(lhs, rhs) result(hyd_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "min_dble_hyd" :: min_dble_hyd

          
          real(PR), intent(in)        :: lhs
          TYPE(hyperdual), intent(in) :: rhs
          TYPE(hyperdual)             :: hyd_out

          if (lhs.LE.rhs) then
            hyd_out = lhs
          else
            hyd_out = rhs
          endif
          
          return
        end function
        

      !----- Maxval
        function maxval_hyd_array(hyd_in) result(hyd_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "maxval_hyd_array" :: maxval_hyd_array

          
          TYPE(hyperdual), dimension(:), intent(in) :: hyd_in
          TYPE(hyperdual)                           :: hyd_out
          integer :: i

          hyd_out = hyd_in(1)
          do i = 1, size(hyd_in)
            hyd_out = max(hyd_out, hyd_in(i))
          enddo

          return
        end function


      !----- Matmul
        function matmul_hyd_array_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "matmul_hyd_array_hyd_matrix" :: matmul_hyd_array_hyd_matrix
        
          TYPE(hyperdual), dimension(:), intent(in)   :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(lhs))     :: res

          res%f0 = matmul(lhs%f0, rhs%f0)
          res%f1 = matmul(lhs%f0, rhs%f1) + matmul(lhs%f1, rhs%f0)
          res%f2 = matmul(lhs%f0, rhs%f2) + matmul(lhs%f2, rhs%f0)
          res%f12 = matmul(lhs%f0, rhs%f12) + matmul(lhs%f1, rhs%f2) + matmul(lhs%f2, rhs%f1) + matmul(lhs%f12, rhs%f0)
            
          return
        end function matmul_hyd_array_hyd_matrix

        function matmul_hyd_array_dble_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "matmul_hyd_array_dble_matrix" :: matmul_hyd_array_dble_matrix
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          real(PR), dimension(:,:), intent(in)      :: rhs
          TYPE(hyperdual), dimension(size(lhs))   :: res

          res%f0 = matmul(lhs%f0, rhs) 
          res%f1 = matmul(lhs%f1, rhs)
          res%f2 = matmul(lhs%f2, rhs)
          res%f12 = matmul(lhs%f12, rhs)

          return
        end function matmul_hyd_array_dble_matrix

        function matmul_hyd_matrix_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "matmul_hyd_matrix_hyd_array" :: matmul_hyd_matrix_hyd_array
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          TYPE(hyperdual), dimension(:), intent(in)   :: rhs
          TYPE(hyperdual), dimension(size(lhs,1))   :: res

          res%f0 = matmul(lhs%f0, rhs%f0)
          res%f1 = matmul(lhs%f0, rhs%f1) + matmul(lhs%f1, rhs%f0)
          res%f2 = matmul(lhs%f0, rhs%f2) + matmul(lhs%f2, rhs%f0)
          res%f12 = matmul(lhs%f0, rhs%f12) + matmul(lhs%f1, rhs%f2) + matmul(lhs%f2, rhs%f1) + matmul(lhs%f12, rhs%f0)

          return
        end function matmul_hyd_matrix_hyd_array 

        function matmul_hyd_matrix_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "matmul_hyd_matrix_hyd_matrix" :: matmul_hyd_matrix_hyd_matrix
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(lhs,1), size(rhs,2)) :: res

          res%f0 = matmul(lhs%f0, rhs%f0)
          res%f1 = matmul(lhs%f0, rhs%f1) + matmul(lhs%f1, rhs%f0)
          res%f2 = matmul(lhs%f0, rhs%f2) + matmul(lhs%f2, rhs%f0)
          res%f12 = matmul(lhs%f0, rhs%f12) + matmul(lhs%f1, rhs%f2) + matmul(lhs%f2, rhs%f1) + matmul(lhs%f12, rhs%f0)

          return

        end function

        function matmul_hyd_matrix_dble_array(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "matmul_hyd_matrix_dble_array" :: matmul_hyd_matrix_dble_array
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          real(PR), dimension(:), intent(in)          :: rhs
          TYPE(hyperdual), dimension(size(lhs,1))   :: res

          res%f0 = matmul(lhs%f0, rhs) 
          res%f1 = matmul(lhs%f1, rhs)
          res%f2 = matmul(lhs%f2, rhs)
          res%f12 = matmul(lhs%f12, rhs)

          return
        end function matmul_hyd_matrix_dble_array 

        function matmul_hyd_matrix_dble_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "matmul_hyd_matrix_dble_matrix" :: matmul_hyd_matrix_dble_matrix
        
          
          TYPE(hyperdual), dimension(:,:), intent(in) :: lhs
          real(PR), dimension(:,:), intent(in)        :: rhs
          TYPE(hyperdual), dimension(size(lhs,1), size(lhs,2))  :: res

          res%f0 = matmul(lhs%f0, rhs)
          res%f1 = matmul(lhs%f1, rhs)
          res%f2 = matmul(lhs%f2, rhs)
          res%f12 = matmul(lhs%f12, rhs)

          return

        end function matmul_hyd_matrix_dble_matrix 

        function matmul_dble_array_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "matmul_dble_array_hyd_matrix" :: matmul_dble_array_hyd_matrix
        
          
          real(PR), dimension(:), intent(in)          :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(lhs))     :: res

          res%f0 = matmul(lhs, rhs%f0)
          res%f1 = matmul(lhs, rhs%f1)
          res%f2 = matmul(lhs, rhs%f2)
          res%f12 = matmul(lhs, rhs%f12)

          return

        end function matmul_dble_array_hyd_matrix 

        function matmul_dble_matrix_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "matmul_dble_matrix_hyd_array" :: matmul_dble_matrix_hyd_array
        
          
          real(PR), dimension(:,:), intent(in)      :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs))  :: res

          res%f0 = matmul(lhs, rhs%f0)
          res%f1 = matmul(lhs, rhs%f1)
          res%f2 = matmul(lhs, rhs%f2)
          res%f12 = matmul(lhs, rhs%f12)

          return

        end function matmul_dble_matrix_hyd_array 

        function matmul_dble_matrix_hyd_matrix(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "matmul_dble_matrix_hyd_matrix" :: matmul_dble_matrix_hyd_matrix
        
          
          real(PR), dimension(:,:), intent(in)        :: lhs
          TYPE(hyperdual), dimension(:,:), intent(in) :: rhs
          TYPE(hyperdual), dimension(size(rhs,1),size(rhs,2))  :: res

          res%f0 = matmul(lhs, rhs%f0)
          res%f1 = matmul(lhs, rhs%f1)
          res%f2 = matmul(lhs, rhs%f2)
          res%f12 = matmul(lhs, rhs%f12)
 
          return

        end function matmul_dble_matrix_hyd_matrix 

      !----- Dot Product
        function dot_product_hyd_array_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "dot_product_hyd_array_hyd_array" :: dot_product_hyd_array_hyd_array
        
          
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual)                           :: res

          res%f0 = dot_product(lhs%f0, rhs%f0)
          res%f1 = dot_product(lhs%f0, rhs%f1) + dot_product(lhs%f1, rhs%f0)
          res%f2 = dot_product(lhs%f0, rhs%f2) + dot_product(lhs%f2, rhs%f0)
          res%f12 = dot_product(lhs%f0, rhs%f12) + dot_product(lhs%f1, rhs%f2) + dot_product(lhs%f2, rhs%f1) &
           + dot_product(lhs%f12, rhs%f0)

          return


        end function dot_product_hyd_array_hyd_array

        function dot_product_hyd_array_dble_array(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "dot_product_hyd_array_dble_array" :: dot_product_hyd_array_dble_array
        
            
          TYPE(hyperdual), dimension(:), intent(in) :: lhs
          real(PR), dimension(:), intent(in)        :: rhs
          TYPE(hyperdual)                           :: res

          res%f0 = dot_product(lhs%f0, rhs) 
          res%f1 = dot_product(lhs%f1, rhs)
          res%f2 = dot_product(lhs%f2, rhs)
          res%f12 = dot_product(lhs%f12, rhs)

          return

        end function dot_product_hyd_array_dble_array
        
        function dot_product_dble_array_hyd_array(lhs, rhs) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "dot_product_dble_array_hyd_array" :: dot_product_dble_array_hyd_array
        
          
          real(PR), dimension(:), intent(in)        :: lhs
          TYPE(hyperdual), dimension(:), intent(in) :: rhs
          TYPE(hyperdual)                           :: res

          res%f0 = dot_product(lhs, rhs%f0)
          res%f1 = dot_product(lhs, rhs%f0)
          res%f2 = dot_product(lhs, rhs%f0)
          res%f12 = dot_product(lhs, rhs%f0)

          return

        end function dot_product_dble_array_hyd_array


      !----- Conjugate
        function conjg_hyd(hyd) result(hydconj)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "conjg_hyd" :: conjg_hyd
          
          
          TYPE(hyperdual), intent(in) :: hyd
          TYPE(hyperdual)             :: hydconj
        
          hydconj%f0 = hyd%f0 
          hydconj%f1 = -hyd%f1
            
          return
        end function conjg_hyd
        
      !----- Real part
        function real_hyd(hyd) result(res)
        
          
          TYPE(hyperdual), intent(in) :: hyd
          real(PR)                    :: res
        
          res = dreal(hyd%f0)
        
          return
        end function real_hyd
        
      !----- Imaginary part
        function imag1_hyd(hyd,iright) result(imag_hyd)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "imag1_hyd" :: imag1_hyd
        
          
          TYPE(hyperdual), intent(in) :: hyd
          integer, intent(in)         :: iright
          real(PR)                    :: imag_hyd
        
          select case (iright)
            case (1)
              imag_hyd  = dimag(hyd%f0)
            case (2)
              imag_hyd  = dreal(hyd%f1)
            case (3)
              imag_hyd  = dimag(hyd%f1)
          end select
        
          return
        end function imag1_hyd
        
        function imag2_array_hyd(hyd,iright) result(qimag)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "imag2_array_hyd" :: imag2_array_hyd
        
          
          TYPE(hyperdual), dimension (:), intent(in) :: hyd
          integer, intent(in)           :: iright
          real(PR), dimension(size(hyd))  :: qimag
        
          select case (iright)
            case (1)
              hydimag = dimag(hyd%f0)
            case (2)
              hydimag = dreal(hyd%f1)
            case (3)
              hydimag = dimag(hyd%f1)
          end select
        
          return
        end function imag2_array_hyd
        
      !----- Exponential
        function exp_hyd(hyd) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "exp_hyd" :: exp_hyd
        
          
          type(hyperdual), intent(in) :: hyd
           type(hyperdual) :: res
           double precision :: deriv
           double precision :: funval
           funval = cos(hyd%f0)
           deriv = -1 * sin(hyd%f0)
           res%f0 = funval
           res%f1 = deriv * hyd%f1
           res%f2 = deriv * hyd%f2
           res%f12 = deriv * hyd%f12 - funval * hyd%f1 * hyd%f2
           return
        end function exp_hyd
        
      !----- Cosine
        function cos_hyd(hyd) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "cos_hyd" :: cos_hyd
        
          
          type(hyperdual), intent(in) :: hyd
          type(hyperdual) :: res
          double precision :: deriv
          double precision :: funval
          funval = cos(hyd%f0)
          deriv = -1 * sin(hyd%f0)
          res%f0 = funval
          res%f1 = deriv * hyd%f1
          res%f2 = deriv * hyd%f2
          res%f12 = deriv * hyd%f12 - funval * hyd%f1 * hyd%f2

          return
            
        end function cos_hyd
      !----- Arccosine
       function acos_hyd(hyd) result(res)
           type(hyperdual), intent(in) :: hyd
           type(hyperdual) :: res
           double precision :: deriv
           double precision :: deriv1
           double precision :: funval
           funval = acos(hyd%f0)
           deriv1 = 1.0 - hyd%f0 * hyd%f0
           deriv = -1.0 / sqrt(deriv1)
           res%f0 = funval
           res%f1 = deriv * hyd%f1
           res%f2 = deriv * hyd%f2
           res%f12 = deriv * hyd%f12 + hyd%f1 * hyd%f2 * ((-1 * hyd%f0)* (deriv1**(-1.5)))
           return
        end function acos_hyd

        
      !----- Sine
        function sin_hyd(hyd) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "sin_hyd" :: sin_hyd
        
          
          type(hyperdual), intent(in) :: hyd
          type(hyperdual) :: res
          double precision :: deriv
          double precision :: funval
          funval = sin(hyd%f0)
          deriv = cos(hyd%f0)
          res%f0 = funval
          res%f1 = deriv * hyd%f1
          res%f2 = deriv * hyd%f2
          res%f12 = deriv * hyd%f12 - funval * hyd%f1 * hyd%f2
          return
          
        end function sin_hyd


      
       
       !----- Arcsine
        function asin_hyd(hyd) result(res)
            type(hyperdual), intent(in) :: hyd
            type(hyperdual) :: res
            double precision :: deriv
            double precision :: deriv1
            double precision :: funval
            funval = asin(hyd%f0)
            deriv1 = 1.0 - hyd%f0 * hyd%f0
            deriv = 1.0/sqrt(deriv1)
            res%f0 = funval
            res%f1 = deriv * hyd%f1
            res%f2 = deriv * hyd%f2
            res%f12 = deriv * hyd%f12 + hyd%f1 * hyd%f2 * (hyd%f0* (deriv1**(-1.5)))
            return
        end function asin_hyd
      !----- Tangent
        function tan_hyd(hyd) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "tan_hyd" :: tan_hyd
        
          
          type(hyperdual), intent(in) :: hyd
          type(hyperdual) :: res
          double precision :: deriv
          double precision :: funval
          funval = tan(hyd%f0)
          deriv = funval * funval + 1.0
          res%f0 = funval
          res%f1 = deriv * hyd%f1
          res%f2 = deriv * hyd%f2
          res%f12 = deriv * hyd%f12 + hyd%f1 * hyd%f2 * (2 * funval * deriv)
          return
        end function tan_hyd

      !----- Arctangent
        function atan_hyd(hyd) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "atan_hyd" :: atan_hyd
          
          type(hyperdual), intent(in) :: hyd
          type(hyperdual) :: res
          double precision :: deriv
          double precision :: deriv1
          double precision :: funval
          funval = atan(hyd%f0)
          deriv1 = 1.0 + hyd%f0 * hyd%f0
          deriv = 1.0 / sqrt(deriv1)
          res%f0 = funval
          res%f1 = deriv * hyd%f1
          res%f2 = deriv * hyd%f2
          res%f12 = deriv * hyd%f12 + hyd%f1 * hyd%f2 * (-2 * hyd%f0 / (deriv1 * deriv1))
          return  
        end function atan_hyd
                
      !Hyperbolic Functions
     
      !Hyperbolic Sine
       function hyd_sinh(hyd) result(res)
            TYPE(hyperdual), intent(in) :: hyd
            TYPE(hyperdual)             :: res
            TYPE(hyperdual)             :: minus
            TYPE(hyperdual)             :: exp1, exp2 
            exp1 = hyd_exp(hyd)
            exp2 = hyd_exp(-hyd)
            minus = hyd_minus_hyd(exp1, exp2)
            res = hyd_div_dble(minus, 2.0)
            
            return
        end function hyd_sinh 
        
      !----- hyperdual module
        function mod_hyd(q) result(r)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "mod_hyd" :: mod_hyd
        
          
          TYPE(hyperdual), intent(in) :: q
          complex(PR)                 :: r
        
          r = sqrt(q%f0**2 + q%f1**2)
        
        end function mod_hyd
        
      !----- hyperdual argument
        function arg_hyd(q) result(theta)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "arg_hyd" :: arg_hyd
        
          
          TYPE(hyperdual), intent(in) :: q
          complex(PR)                 :: theta
        
            theta = cdatan2(q%f1, q%f0)
            
        end function arg_hyd
        
      !----- hyperdual logarithm
        function log_hyd(hyd) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "log_hyd" :: log_hyd
           type(hyperdual), intent(in) :: hyd
           type(hyperdual) :: res
           double precision :: deriv1
           double precision :: deriv2
           deriv1 = hyd%f1/hyd%f0
           deriv2 = hyd%f2/hyd%f0
           res%f0 = log(hyd%f0)
           res%f1 = deriv1
           res%f2 = deriv2
           res%f12 = hyd%f12/hyd%f0 - (deriv1 * deriv2)
           return
        end function log_hyd
        
      !----- hyperdual square root
        function sqrt_hyd(hyd) result(res)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "sqrt_hyd" :: sqrt_hyd
           type(hyperdual), intent(in) :: hyd
           type(hyperdual) :: res
           res = hyd_pow_dble(hyd, 0.5)
        
        end function sqrt_hyd
        
      !----- hyperdual log
        !function bdlog1p(z) result(zlog1p)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "bdlog1p" :: bdlog1p
            
         ! 
          
           !TYPE(hyperdual), intent(in) :: z
          !TYPE(hyperdual)             :: zlog1p
          !complex(PR)                 :: DATAN2
          !complex(PR)                 :: x,y
          !real(PR)                    :: d1mach
        
        !  x = z%f0
         ! y = z%f1
          !if (abs(real(x,PR))*d1mach(4) <= 1.0_PR) then
           ! zlog1p = hyd(0.5_PR*CDLOG1P(x * (2.0_PR + x) + y*y),CDATAN2(y,1.0_PR + x))
          !else
          !  zlog1p = BDLOG(z)
          !end if
            
        !end function
        

      !----------------------------------------------------------------!
      !                     COMPARISON OPERATORS                       !
      !----------------------------------------------------------------!

      

      !----------------------------------------------------------------!
      !                     INTRINSIC FUNCTIONS                        !
      !----------------------------------------------------------------!

      !----- ABS (absolute value)
        function abs_cplx(hyd_in) result(X_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "abs_cplx" :: abs_cplx
          
          
          complex(PR), intent(in)   :: hyd_in
          complex(PR)               :: X_out

          X_out = sign(1.d0, REAL(hyd_in,PR)) * hyd_in 
        
        end function

        function abs_cplx_array(hyd_in) result(X_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "abs_cplx_array" :: abs_cplx_array
          
          
          complex(PR), dimension(:), intent(in) :: hyd_in
          complex(PR), dimension(size(hyd_in))    :: X_out
          integer :: i
    
          do i = 1, size(hyd_in)
            X_out(i) = abs(hyd_in(i))
          enddo

        end function

        function abs_cplx_matrix(hyd_in) result(X_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "abs_cplx_matrix" :: abs_cplx_matrix
          
          
          complex(PR), dimension(:,:), intent(in)           :: hyd_in
          complex(PR), dimension(size(hyd_in,1), size(hyd_in,2)):: X_out
          integer :: i, j

          do i = 1, size(hyd_in,1)
            do j = 1, size(hyd_in,2)
              X_out(i,j) = abs(hyd_in(i,j))
            enddo
          enddo

        end function


      !----- SIGN
        function sign_cplx_cplx(lhs, rhs) result(hyd_out)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : "sign_cplx_cplx" :: sign_cplx_cplx
          
          
          complex(PR), intent(in) :: lhs, rhs
          complex(PR)             :: hyd_out

          if (REAL(rhs,PR).GE.(0.d0)) then
            hyd_out = abs(lhs)
          else
            hyd_out = -abs(lhs)
          endif

        end function sign_cplx_cplx
    end module hydfinal