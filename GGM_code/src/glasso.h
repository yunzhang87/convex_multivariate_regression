#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif
    
    
    enum penalty_types
    {
        L1,
        Truncated_L1,
        MCP
    };
    
    
    typedef struct
    {
        double  *Delta, *Delta_old, *tmp1, *tmp2, *Lambda, *lam_mat;
    } tmpvars;
    
    
    typedef struct
    {
        double  *Delta, *Delta_old, *tmp1, *tmp2, *Lambda, *lam_mat_L1, *lam_mat_group;
    } tmpvars_group;
    
    
    
    
#ifdef __cplusplus
}
#endif
