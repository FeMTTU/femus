set(rifat_SRCS
    rifat.cpp
)

ecm_qt_declare_logging_category(rifat_SRCS
    HEADER debug.h
    IDENTIFIER PLUGIN_RIFAT
    CATEGORY_NAME "kdevelop.plugins.rifat"
)

kdevplatform_add_plugin(rifat
    JSON rifat.json
    SOURCES ${rifat_SRCS}
)

target_link_libraries(rifat
    KDev::Interfaces
)
