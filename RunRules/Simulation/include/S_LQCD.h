/**\file */
#ifndef SLIC_DECLARATIONS_S_LQCD_H
#define SLIC_DECLARATIONS_S_LQCD_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define S_LQCD_LZ (16)
#define S_LQCD_LY (16)
#define S_LQCD_LX (32)
#define S_LQCD_T (128)
#define S_LQCD_gCmdSize (128)
#define S_LQCD_spCmdSize (128)
#define S_LQCD_loopOffset (16)
#define S_LQCD_numPipes (8)
#define S_LQCD_PCIE_ALIGNMENT (16)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] ticks_diracKernel The number of ticks for which kernel "diracKernel" will run.
 * \param [in] ticks_gReadCmdKernel The number of ticks for which kernel "gReadCmdKernel" will run.
 * \param [in] ticks_gWriteCmdKernel The number of ticks for which kernel "gWriteCmdKernel" will run.
 * \param [in] ticks_spReadCmdKernel0 The number of ticks for which kernel "spReadCmdKernel0" will run.
 * \param [in] ticks_spReadCmdKernel1 The number of ticks for which kernel "spReadCmdKernel1" will run.
 * \param [in] ticks_spWriteCmdKernel The number of ticks for which kernel "spWriteCmdKernel" will run.
 * \param [in] inscalar_diracKernel_alpha Input scalar parameter "diracKernel.alpha".
 * \param [in] inscalar_diracKernel_beta_s Input scalar parameter "diracKernel.beta_s".
 * \param [in] inscalar_diracKernel_beta_t_b Input scalar parameter "diracKernel.beta_t_b".
 * \param [in] inscalar_diracKernel_beta_t_f Input scalar parameter "diracKernel.beta_t_f".
 * \param [in] inscalar_diracKernel_doSub Input scalar parameter "diracKernel.doSub".
 * \param [in] inscalar_diracKernel_ieo Input scalar parameter "diracKernel.ieo".
 * \param [in] inscalar_diracKernel_isign Input scalar parameter "diracKernel.isign".
 * \param [in] inscalar_gReadCmdKernel_startAddress Input scalar parameter "gReadCmdKernel.startAddress".
 * \param [in] inscalar_gWriteCmdKernel_startAddress Input scalar parameter "gWriteCmdKernel.startAddress".
 * \param [in] inscalar_spReadCmdKernel0_halos Input scalar parameter "spReadCmdKernel0.halos".
 * \param [in] inscalar_spReadCmdKernel0_startAddress Input scalar parameter "spReadCmdKernel0.startAddress".
 * \param [in] inscalar_spReadCmdKernel1_halos Input scalar parameter "spReadCmdKernel1.halos".
 * \param [in] inscalar_spReadCmdKernel1_startAddress Input scalar parameter "spReadCmdKernel1.startAddress".
 * \param [in] inscalar_spWriteCmdKernel_halos Input scalar parameter "spWriteCmdKernel.halos".
 * \param [in] inscalar_spWriteCmdKernel_startAddress Input scalar parameter "spWriteCmdKernel.startAddress".
 * \param [in] instream_gauge_in Stream "gauge_in".
 * \param [in] instream_size_gauge_in The size of the stream instream_gauge_in in bytes.
 * \param [in] instream_spinor_in Stream "spinor_in".
 * \param [in] instream_size_spinor_in The size of the stream instream_spinor_in in bytes.
 * \param [out] outstream_spinor_out Stream "spinor_out".
 * \param [in] outstream_size_spinor_out The size of the stream outstream_spinor_out in bytes.
 * \param [in] routing_string A string containing comma-separated "from_name -> to_name" routing commands.
 */
void S_LQCD(
	uint64_t ticks_diracKernel,
	uint64_t ticks_gReadCmdKernel,
	uint64_t ticks_gWriteCmdKernel,
	uint64_t ticks_spReadCmdKernel0,
	uint64_t ticks_spReadCmdKernel1,
	uint64_t ticks_spWriteCmdKernel,
	double inscalar_diracKernel_alpha,
	double inscalar_diracKernel_beta_s,
	double inscalar_diracKernel_beta_t_b,
	double inscalar_diracKernel_beta_t_f,
	uint64_t inscalar_diracKernel_doSub,
	uint64_t inscalar_diracKernel_ieo,
	uint64_t inscalar_diracKernel_isign,
	uint64_t inscalar_gReadCmdKernel_startAddress,
	uint64_t inscalar_gWriteCmdKernel_startAddress,
	uint64_t inscalar_spReadCmdKernel0_halos,
	uint64_t inscalar_spReadCmdKernel0_startAddress,
	uint64_t inscalar_spReadCmdKernel1_halos,
	uint64_t inscalar_spReadCmdKernel1_startAddress,
	uint64_t inscalar_spWriteCmdKernel_halos,
	uint64_t inscalar_spWriteCmdKernel_startAddress,
	const void *instream_gauge_in,
	size_t instream_size_gauge_in,
	const void *instream_spinor_in,
	size_t instream_size_spinor_in,
	void *outstream_spinor_out,
	size_t outstream_size_spinor_out,
	const char * routing_string);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] ticks_diracKernel The number of ticks for which kernel "diracKernel" will run.
 * \param [in] ticks_gReadCmdKernel The number of ticks for which kernel "gReadCmdKernel" will run.
 * \param [in] ticks_gWriteCmdKernel The number of ticks for which kernel "gWriteCmdKernel" will run.
 * \param [in] ticks_spReadCmdKernel0 The number of ticks for which kernel "spReadCmdKernel0" will run.
 * \param [in] ticks_spReadCmdKernel1 The number of ticks for which kernel "spReadCmdKernel1" will run.
 * \param [in] ticks_spWriteCmdKernel The number of ticks for which kernel "spWriteCmdKernel" will run.
 * \param [in] inscalar_diracKernel_alpha Input scalar parameter "diracKernel.alpha".
 * \param [in] inscalar_diracKernel_beta_s Input scalar parameter "diracKernel.beta_s".
 * \param [in] inscalar_diracKernel_beta_t_b Input scalar parameter "diracKernel.beta_t_b".
 * \param [in] inscalar_diracKernel_beta_t_f Input scalar parameter "diracKernel.beta_t_f".
 * \param [in] inscalar_diracKernel_doSub Input scalar parameter "diracKernel.doSub".
 * \param [in] inscalar_diracKernel_ieo Input scalar parameter "diracKernel.ieo".
 * \param [in] inscalar_diracKernel_isign Input scalar parameter "diracKernel.isign".
 * \param [in] inscalar_gReadCmdKernel_startAddress Input scalar parameter "gReadCmdKernel.startAddress".
 * \param [in] inscalar_gWriteCmdKernel_startAddress Input scalar parameter "gWriteCmdKernel.startAddress".
 * \param [in] inscalar_spReadCmdKernel0_halos Input scalar parameter "spReadCmdKernel0.halos".
 * \param [in] inscalar_spReadCmdKernel0_startAddress Input scalar parameter "spReadCmdKernel0.startAddress".
 * \param [in] inscalar_spReadCmdKernel1_halos Input scalar parameter "spReadCmdKernel1.halos".
 * \param [in] inscalar_spReadCmdKernel1_startAddress Input scalar parameter "spReadCmdKernel1.startAddress".
 * \param [in] inscalar_spWriteCmdKernel_halos Input scalar parameter "spWriteCmdKernel.halos".
 * \param [in] inscalar_spWriteCmdKernel_startAddress Input scalar parameter "spWriteCmdKernel.startAddress".
 * \param [in] instream_gauge_in Stream "gauge_in".
 * \param [in] instream_size_gauge_in The size of the stream instream_gauge_in in bytes.
 * \param [in] instream_spinor_in Stream "spinor_in".
 * \param [in] instream_size_spinor_in The size of the stream instream_spinor_in in bytes.
 * \param [out] outstream_spinor_out Stream "spinor_out".
 * \param [in] outstream_size_spinor_out The size of the stream outstream_spinor_out in bytes.
 * \param [in] routing_string A string containing comma-separated "from_name -> to_name" routing commands.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *S_LQCD_nonblock(
	uint64_t ticks_diracKernel,
	uint64_t ticks_gReadCmdKernel,
	uint64_t ticks_gWriteCmdKernel,
	uint64_t ticks_spReadCmdKernel0,
	uint64_t ticks_spReadCmdKernel1,
	uint64_t ticks_spWriteCmdKernel,
	double inscalar_diracKernel_alpha,
	double inscalar_diracKernel_beta_s,
	double inscalar_diracKernel_beta_t_b,
	double inscalar_diracKernel_beta_t_f,
	uint64_t inscalar_diracKernel_doSub,
	uint64_t inscalar_diracKernel_ieo,
	uint64_t inscalar_diracKernel_isign,
	uint64_t inscalar_gReadCmdKernel_startAddress,
	uint64_t inscalar_gWriteCmdKernel_startAddress,
	uint64_t inscalar_spReadCmdKernel0_halos,
	uint64_t inscalar_spReadCmdKernel0_startAddress,
	uint64_t inscalar_spReadCmdKernel1_halos,
	uint64_t inscalar_spReadCmdKernel1_startAddress,
	uint64_t inscalar_spWriteCmdKernel_halos,
	uint64_t inscalar_spWriteCmdKernel_startAddress,
	const void *instream_gauge_in,
	size_t instream_size_gauge_in,
	const void *instream_spinor_in,
	size_t instream_size_spinor_in,
	void *outstream_spinor_out,
	size_t outstream_size_spinor_out,
	const char * routing_string);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint64_t ticks_diracKernel; /**<  [in] The number of ticks for which kernel "diracKernel" will run. */
	uint64_t ticks_gReadCmdKernel; /**<  [in] The number of ticks for which kernel "gReadCmdKernel" will run. */
	uint64_t ticks_gWriteCmdKernel; /**<  [in] The number of ticks for which kernel "gWriteCmdKernel" will run. */
	uint64_t ticks_spReadCmdKernel0; /**<  [in] The number of ticks for which kernel "spReadCmdKernel0" will run. */
	uint64_t ticks_spReadCmdKernel1; /**<  [in] The number of ticks for which kernel "spReadCmdKernel1" will run. */
	uint64_t ticks_spWriteCmdKernel; /**<  [in] The number of ticks for which kernel "spWriteCmdKernel" will run. */
	double inscalar_diracKernel_alpha; /**<  [in] Input scalar parameter "diracKernel.alpha". */
	double inscalar_diracKernel_beta_s; /**<  [in] Input scalar parameter "diracKernel.beta_s". */
	double inscalar_diracKernel_beta_t_b; /**<  [in] Input scalar parameter "diracKernel.beta_t_b". */
	double inscalar_diracKernel_beta_t_f; /**<  [in] Input scalar parameter "diracKernel.beta_t_f". */
	uint64_t inscalar_diracKernel_doSub; /**<  [in] Input scalar parameter "diracKernel.doSub". */
	uint64_t inscalar_diracKernel_ieo; /**<  [in] Input scalar parameter "diracKernel.ieo". */
	uint64_t inscalar_diracKernel_isign; /**<  [in] Input scalar parameter "diracKernel.isign". */
	uint64_t inscalar_gReadCmdKernel_startAddress; /**<  [in] Input scalar parameter "gReadCmdKernel.startAddress". */
	uint64_t inscalar_gWriteCmdKernel_startAddress; /**<  [in] Input scalar parameter "gWriteCmdKernel.startAddress". */
	uint64_t inscalar_spReadCmdKernel0_halos; /**<  [in] Input scalar parameter "spReadCmdKernel0.halos". */
	uint64_t inscalar_spReadCmdKernel0_startAddress; /**<  [in] Input scalar parameter "spReadCmdKernel0.startAddress". */
	uint64_t inscalar_spReadCmdKernel1_halos; /**<  [in] Input scalar parameter "spReadCmdKernel1.halos". */
	uint64_t inscalar_spReadCmdKernel1_startAddress; /**<  [in] Input scalar parameter "spReadCmdKernel1.startAddress". */
	uint64_t inscalar_spWriteCmdKernel_halos; /**<  [in] Input scalar parameter "spWriteCmdKernel.halos". */
	uint64_t inscalar_spWriteCmdKernel_startAddress; /**<  [in] Input scalar parameter "spWriteCmdKernel.startAddress". */
	const void *instream_gauge_in; /**<  [in] Stream "gauge_in". */
	size_t instream_size_gauge_in; /**<  [in] The size of the stream instream_gauge_in in bytes. */
	const void *instream_spinor_in; /**<  [in] Stream "spinor_in". */
	size_t instream_size_spinor_in; /**<  [in] The size of the stream instream_spinor_in in bytes. */
	void *outstream_spinor_out; /**<  [out] Stream "spinor_out". */
	size_t outstream_size_spinor_out; /**<  [in] The size of the stream outstream_spinor_out in bytes. */
	const char * routing_string; /**<  [in] A string containing comma-separated "from_name -> to_name" routing commands. */
} S_LQCD_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void S_LQCD_run(
	max_engine_t *engine,
	S_LQCD_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *S_LQCD_run_nonblock(
	max_engine_t *engine,
	S_LQCD_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void S_LQCD_run_group(max_group_t *group, S_LQCD_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *S_LQCD_run_group_nonblock(max_group_t *group, S_LQCD_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void S_LQCD_run_array(max_engarray_t *engarray, S_LQCD_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *S_LQCD_run_array_nonblock(max_engarray_t *engarray, S_LQCD_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* S_LQCD_convert(max_file_t *maxfile, S_LQCD_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* S_LQCD_init(void);

/* Error handling functions */
int S_LQCD_has_errors(void);
const char* S_LQCD_get_errors(void);
void S_LQCD_clear_errors(void);
/* Free statically allocated maxfile data */
void S_LQCD_free(void);
/* returns: -1 = error running command; 0 = no error reported */
int S_LQCD_simulator_start(void);
/* returns: -1 = error running command; 0 = no error reported */
int S_LQCD_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_S_LQCD_H */

