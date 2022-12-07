/**
 * @file new_geometry.h
 * @brief Functions needed for the new geometry implementation that
 *        will replace the current geometry in the future
 */

 #ifndef NEW_GEOMETRY_H
 #define NEW_GEOMETRY_H
 #ifdef __cplusplus
    extern "C" {
 #endif

void define_geometry();
void* sendbuf_alloc(size_t bytes_per_site);
void sync_field(geometry_descriptor *gd, int byte_per_site, int is_spinor_like, void *latticebuf, void *sb_ptr);
int test_define_geometry();
void sendbuf_report();

#ifdef __cplusplus
    }
#endif
 #endif