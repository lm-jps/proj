c#######################################################################
      module syntax
c
c-----------------------------------------------------------------------
c ****** Group definitions for parsing command-line arguments.
c-----------------------------------------------------------------------
c
c        GROUP 1: <kw>
c        GROUP 2: <arg>
c        GROUP 3: <kw> <arg>
c        GROUP 4: <kw> <arg> <arg>
c
      integer, parameter :: ngroups=4
c
      integer, parameter :: GROUP_K  =1
      integer, parameter :: GROUP_A  =2
      integer, parameter :: GROUP_KA =3
      integer, parameter :: GROUP_KAA=4
c
      end module
c#######################################################################
      module string_def
c
c-----------------------------------------------------------------------
c ****** Define a structure to hold a string.
c-----------------------------------------------------------------------
c
      implicit none
c
      type :: string
        character, dimension(:), pointer :: c
      end type
c
      end module
c#######################################################################
      module paragraph_def
c
      use string_def
c
      implicit none
c
c-----------------------------------------------------------------------
c ****** Define a structure for a linked list of lines
c ****** (i.e., a paragraph).
c-----------------------------------------------------------------------
c
      type :: paragraph
        type(string) :: line
        type(paragraph), pointer :: next
      end type
c
c-----------------------------------------------------------------------
c ****** Define a structure to hold a list of paragraphs.
c-----------------------------------------------------------------------
c
      type :: parlist
        type(paragraph), pointer :: par
      end type
c
      end module
c#######################################################################
      module lcase_interface
      interface
        function lcase (s)
        character(*) :: s
        character(len(s)) :: lcase
        end
      end interface
      end module
c#######################################################################
      module ucase_interface
      interface
        function ucase (s)
        character(*) :: s
        character(len(s)) :: ucase
        end
      end interface
      end module
c#######################################################################
      module new_par_interface
      interface
        subroutine new_par (par)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: par
        end
      end interface
      end module
c#######################################################################
      module delete_par_interface
      interface
        subroutine delete_par (par)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: par
        end
      end interface
      end module
c#######################################################################
      module add_line_interface
      interface
        subroutine add_line (line,par)
        use paragraph_def
        implicit none
        character(*) :: line
        type(paragraph), pointer :: par
        end
      end interface
      end module
c#######################################################################
      module print_par_interface
      interface
        subroutine print_par (par)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: par
        end
      end interface
      end module
c#######################################################################
      module get_str_interface
      interface
        function get_str (str)
        use string_def
        implicit none
        type(string) :: str
        character(size(str%c)) :: get_str
        end
      end interface
      end module
c#######################################################################
      module get_usage_line_interface
      interface
        subroutine get_usage_line (usage)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: usage
        end
      end interface
      end module
