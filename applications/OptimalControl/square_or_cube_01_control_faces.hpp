#ifndef CONTROL_FACES_SQUARE_OR_CUBE_HPP
#define CONTROL_FACES_SQUARE_OR_CUBE_HPP


namespace femus {

namespace ctrl {


   class Domain_of_PDE {

  public:


 };
 
 

 namespace square_or_cube {
     

   class List_of_faces : public Domain_of_PDE {

  public:


   static  const  int sign_function_for_delimiting_region(const unsigned int face_index) {

   int  target_line_sign;

        if (face_index == 1 || face_index == 3 || face_index == 5) { target_line_sign = 1;  }
   else if (face_index == 2 || face_index == 4 || face_index == 6) { target_line_sign = -1; }

   return target_line_sign;

   }


 };
 
 //*********************** Gamma_c, list of control faces - BEGIN *****************************************

   class List_of_Gamma_control_faces : public List_of_faces {

  public:



 //direction of the line that contains \Gamma_c
      static const unsigned int normal_direction_to_Gamma_control(const unsigned int face_index) {

    unsigned int axis_dir;

        if (face_index == 1 || face_index == 2) { axis_dir = 0; }
   else if (face_index == 3 || face_index == 4) { axis_dir = 1; }
   else if (face_index == 5 || face_index == 6) { axis_dir = 2; }

    return axis_dir;

}




      static const unsigned int tangential_direction_to_Gamma_control(const unsigned int face_index) {

    unsigned int axis_dir;

        if (face_index == 1 || face_index == 2) { axis_dir = 1; }
   else if (face_index == 3 || face_index == 4) { axis_dir = 0; }
   else if (face_index == 5 || face_index == 6) { /*abort();*/ axis_dir = 1; }  ///@todo this depends on the mesh file

    return axis_dir;

}


      static const unsigned int opposite_face(const unsigned int face_index) {

   unsigned int face_opposite = 0;

        if (face_index == 1 || face_index == 3 || face_index == 5) { face_opposite = face_index + 1; }
   else if (face_index == 2 || face_index == 4 || face_index == 6) { face_opposite = face_index - 1; }

    return face_opposite;

}

   };


//*********************** Single - BEGIN *****************************************


  class List_of_Gamma_control_faces_Single : public List_of_Gamma_control_faces {

  public:

     static constexpr unsigned _face_with_extremes_index_size = 1 ;

     static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
         true
    };

     static const double   _face_with_extremes_extremes[ _face_with_extremes_index_size ][2];


  };

      const double   List_of_Gamma_control_faces_Single::_face_with_extremes_extremes[ _face_with_extremes_index_size ][2] = {
       { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
    };



  class List_of_Gamma_control_faces_One :  public  List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

     const unsigned List_of_Gamma_control_faces_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    1    };



  class List_of_Gamma_control_faces_Two : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    2    };


  class List_of_Gamma_control_faces_Three : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    3    };



     class List_of_Gamma_control_faces_Four : public List_of_Gamma_control_faces_Single {

  public:

     static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ];

  };

      const unsigned List_of_Gamma_control_faces_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Single::_face_with_extremes_index_size ] = {    4    };

//*********************** Single - END *****************************************


//*********************** Double - BEGIN *****************************************


  class List_of_Gamma_control_faces_Double : public List_of_Gamma_control_faces {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 2 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
      };

       static const double   _face_with_extremes_extremes[ _face_with_extremes_index_size ][2];


    };


      const double   List_of_Gamma_control_faces_Double::_face_with_extremes_extremes[ _face_with_extremes_index_size ][2] = {
         { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
       , { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }

      };

      class List_of_Gamma_control_faces_One_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  2   };

    class List_of_Gamma_control_faces_One_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  3   };




     class List_of_Gamma_control_faces_One_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    1  ,  4   };


    class List_of_Gamma_control_faces_Two_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  1   };


    class List_of_Gamma_control_faces_Two_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  3   };



    class List_of_Gamma_control_faces_Two_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    2  ,  4   };




    class List_of_Gamma_control_faces_Three_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  1   };



    class List_of_Gamma_control_faces_Three_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  2   };





    class List_of_Gamma_control_faces_Three_Four :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    3  ,  4   };




    class List_of_Gamma_control_faces_Four_One :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  1   };



    class List_of_Gamma_control_faces_Four_Two :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  2   };





    class List_of_Gamma_control_faces_Four_Three :  public  List_of_Gamma_control_faces_Double {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Double::_face_with_extremes_index_size ] = {    4  ,  3   };




//*********************** Double - END *****************************************


//*********************** Triple - BEGIN *****************************************


  class List_of_Gamma_control_faces_Triple : public List_of_Gamma_control_faces {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 3 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
          ,true
      };

       static const double   _face_with_extremes_extremes[ _face_with_extremes_index_size ][2];


    };


      const double   List_of_Gamma_control_faces_Triple::_face_with_extremes_extremes[ _face_with_extremes_index_size ][2] = {
         { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
       , { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
       , { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }

      };


      class List_of_Gamma_control_faces_One_Three_Two :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  3  ,    2  };



      class List_of_Gamma_control_faces_One_Four_Two :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  4  ,    2  };



      class List_of_Gamma_control_faces_One_Three_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  3  ,    4  };



      class List_of_Gamma_control_faces_One_Four_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    1  ,  4  ,    3  };





      class List_of_Gamma_control_faces_Two_Three_One :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  3  ,    1  };



      class List_of_Gamma_control_faces_Two_Four_One :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  4  ,    1  };



      class List_of_Gamma_control_faces_Two_Three_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Three_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  3  ,    4  };


      class List_of_Gamma_control_faces_Two_Four_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Two_Four_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    2  ,  4  ,    3  };



      class List_of_Gamma_control_faces_Three_One_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  1  ,    4  };



      class List_of_Gamma_control_faces_Three_Two_Four :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two_Four::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  2  ,    4  };


      class List_of_Gamma_control_faces_Three_One_Two:  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  1  ,    2  };


      class List_of_Gamma_control_faces_Three_Two_One:  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Three_Two_One::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    3  ,  2  ,    1  };


      class List_of_Gamma_control_faces_Four_One_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    4  ,  1  ,    3  };



      class List_of_Gamma_control_faces_Four_Two_Three :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_Two_Three::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    4  ,  2  ,    3  };



      class List_of_Gamma_control_faces_Four_One_Two :  public  List_of_Gamma_control_faces_Triple {

        public:

           static const unsigned _face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_Four_One_Two::_face_with_extremes_index[ List_of_Gamma_control_faces_Triple::_face_with_extremes_index_size ] = {    4  ,  1  ,    2  };

//*********************** Triple - END *****************************************


//*********************** Quadruple - BEGIN *****************************************


  class List_of_Gamma_control_faces_Quadruple : public List_of_Gamma_control_faces {

    public:

       static constexpr unsigned _face_with_extremes_index_size = 4 ;

       static constexpr bool     _face_with_extremes_extract_subface[ _face_with_extremes_index_size ] = {
           true
          ,true
          ,true
          ,true
      };

       static const double   _face_with_extremes_extremes[ _face_with_extremes_index_size ][2];


    };


      const double   List_of_Gamma_control_faces_Quadruple::_face_with_extremes_extremes[ _face_with_extremes_index_size ][2] = {
         { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
       , { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
       , { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
       , { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }

      };




      class List_of_Gamma_control_faces_One_Three_Two_Four :  public  List_of_Gamma_control_faces_Quadruple {

        public:

           static const unsigned _face_with_extremes_index[ /*List_of_Gamma_control_faces_Quadruple::*/_face_with_extremes_index_size ];

        };

           const unsigned List_of_Gamma_control_faces_One_Three_Two_Four::_face_with_extremes_index[ /*List_of_Gamma_control_faces_Quadruple::*/_face_with_extremes_index_size ] = {    1  ,  3  , 4  ,   2  };



//*********************** Quadruple - END *****************************************


//*********************** Gamma_c, list of control faces - END *****************************************

 }
 
 
  }

}


#endif

