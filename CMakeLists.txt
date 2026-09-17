cmake_minimum_required(VERSION 3.16)

project(csm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

function(target_copy_bathymetry_data TARGET_NAME)
    set_property(TARGET ${TARGET_NAME} APPEND PROPERTY ASSETS_TO_COPY "data/bathymetry")
endfunction()

function(target_copy_wind_data TARGET_NAME)
    set_property(TARGET ${TARGET_NAME} APPEND PROPERTY ASSETS_TO_COPY "data/wind")
endfunction()

function(target_finalize_asset_copy TARGET_NAME)
    get_target_property(ASSET_LIST ${TARGET_NAME} ASSETS_TO_COPY)
    
    if(ASSET_LIST)
        foreach(ASSET_DIR ${ASSET_LIST})
            list(APPEND COPY_COMMANDS 
                 COMMAND ${CMAKE_COMMAND} -E copy_directory
                         "${CMAKE_SOURCE_DIR}/${ASSET_DIR}"
                         "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ASSET_DIR}"
            )
        endforeach()

        add_custom_command(
            TARGET ${TARGET_NAME}
            POST_BUILD
            ${COPY_COMMANDS}
            COMMENT "Copying registered assets for ${TARGET_NAME}"
        )
    endif()
endfunction()

function(target_add_manifest TARGET_NAME)
    if(MSVC)
        target_sources(${TARGET_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/app.manifest")
    endif()
endfunction()

add_subdirectory(utils)

add_subdirectory(slae)
add_subdirectory(calc)
add_subdirectory(generators)

add_subdirectory(tests)