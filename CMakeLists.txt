option(OSPRAY_MODULE_EXAJET_IMPORTER "Build NASA exajet scenegraph importer module" OFF)


if (OSPRAY_MODULE_EXAJET_IMPORTER)
  ospray_create_library(ospray_module_exajet_import
    import_exajet.cpp
    TAMRLevelKDT.cpp
    3rd_lib/chull.cpp
  LINK
    ospray_sg
    ospray_common
  )
endif()

