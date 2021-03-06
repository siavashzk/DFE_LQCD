/**\file */
#ifndef SLIC_DECLARATIONS_LQCD_H
#define SLIC_DECLARATIONS_LQCD_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define LQCD_LZ (8)
#define LQCD_LY (8)
#define LQCD_LX (16)
#define LQCD_T (8)
#define LQCD_gCmdSize (128)
#define LQCD_spCmdSize (128)
#define LQCD_loopOffset (16)
#define LQCD_numPipes (8)
#define LQCD_PCIE_ALIGNMENT (16)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] ticks_diracKernel The number of ticks for which kernel "diracKernel" will run.
 * \param [in] ticks_readCmdKernel0 The number of ticks for which kernel "readCmdKernel0" will run.
 * \param [in] ticks_readCmdKernel1 The number of ticks for which kernel "readCmdKernel1" will run.
 * \param [in] ticks_readCmdKernel2 The number of ticks for which kernel "readCmdKernel2" will run.
 * \param [in] ticks_writeCmdKernel The number of ticks for which kernel "writeCmdKernel" will run.
 * \param [in] inscalar_diracKernel_alpha Input scalar parameter "diracKernel.alpha".
 * \param [in] inscalar_diracKernel_beta_s Input scalar parameter "diracKernel.beta_s".
 * \param [in] inscalar_diracKernel_beta_t_b Input scalar parameter "diracKernel.beta_t_b".
 * \param [in] inscalar_diracKernel_beta_t_f Input scalar parameter "diracKernel.beta_t_f".
 * \param [in] inscalar_diracKernel_doSub Input scalar parameter "diracKernel.doSub".
 * \param [in] inscalar_diracKernel_ieo Input scalar parameter "diracKernel.ieo".
 * \param [in] inscalar_diracKernel_isign Input scalar parameter "diracKernel.isign".
 * \param [in] inscalar_readCmdKernel0_burstsPerTSlice Input scalar parameter "readCmdKernel0.burstsPerTSlice".
 * \param [in] inscalar_readCmdKernel0_cmdSize Input scalar parameter "readCmdKernel0.cmdSize".
 * \param [in] inscalar_readCmdKernel0_halos Input scalar parameter "readCmdKernel0.halos".
 * \param [in] inscalar_readCmdKernel0_startAddress Input scalar parameter "readCmdKernel0.startAddress".
 * \param [in] inscalar_readCmdKernel1_burstsPerTSlice Input scalar parameter "readCmdKernel1.burstsPerTSlice".
 * \param [in] inscalar_readCmdKernel1_cmdSize Input scalar parameter "readCmdKernel1.cmdSize".
 * \param [in] inscalar_readCmdKernel1_halos Input scalar parameter "readCmdKernel1.halos".
 * \param [in] inscalar_readCmdKernel1_startAddress Input scalar parameter "readCmdKernel1.startAddress".
 * \param [in] inscalar_readCmdKernel2_burstsPerTSlice Input scalar parameter "readCmdKernel2.burstsPerTSlice".
 * \param [in] inscalar_readCmdKernel2_cmdSize Input scalar parameter "readCmdKernel2.cmdSize".
 * \param [in] inscalar_readCmdKernel2_halos Input scalar parameter "readCmdKernel2.halos".
 * \param [in] inscalar_readCmdKernel2_startAddress Input scalar parameter "readCmdKernel2.startAddress".
 * \param [in] inscalar_writeCmdKernel_burstsPerTSlice Input scalar parameter "writeCmdKernel.burstsPerTSlice".
 * \param [in] inscalar_writeCmdKernel_cmdSize Input scalar parameter "writeCmdKernel.cmdSize".
 * \param [in] inscalar_writeCmdKernel_halos Input scalar parameter "writeCmdKernel.halos".
 * \param [in] inscalar_writeCmdKernel_startAddress Input scalar parameter "writeCmdKernel.startAddress".
 * \param [in] instream_data_in Stream "data_in".
 * \param [in] instream_size_data_in The size of the stream instream_data_in in bytes.
 * \param [out] outstream_data_out Stream "data_out".
 * \param [in] outstream_size_data_out The size of the stream outstream_data_out in bytes.
 * \param [in] routing_string A string containing comma-separated "from_name -> to_name" routing commands.
 */
void LQCD(
	uint64_t ticks_diracKernel,
	uint64_t ticks_readCmdKernel0,
	uint64_t ticks_readCmdKernel1,
	uint64_t ticks_readCmdKernel2,
	uint64_t ticks_writeCmdKernel,
	double inscalar_diracKernel_alpha,
	double inscalar_diracKernel_beta_s,
	double inscalar_diracKernel_beta_t_b,
	double inscalar_diracKernel_beta_t_f,
	uint64_t inscalar_diracKernel_doSub,
	uint64_t inscalar_diracKernel_ieo,
	uint64_t inscalar_diracKernel_isign,
	uint64_t inscalar_readCmdKernel0_burstsPerTSlice,
	uint64_t inscalar_readCmdKernel0_cmdSize,
	uint64_t inscalar_readCmdKernel0_halos,
	uint64_t inscalar_readCmdKernel0_startAddress,
	uint64_t inscalar_readCmdKernel1_burstsPerTSlice,
	uint64_t inscalar_readCmdKernel1_cmdSize,
	uint64_t inscalar_readCmdKernel1_halos,
	uint64_t inscalar_readCmdKernel1_startAddress,
	uint64_t inscalar_readCmdKernel2_burstsPerTSlice,
	uint64_t inscalar_readCmdKernel2_cmdSize,
	uint64_t inscalar_readCmdKernel2_halos,
	uint64_t inscalar_readCmdKernel2_startAddress,
	uint64_t inscalar_writeCmdKernel_burstsPerTSlice,
	uint64_t inscalar_writeCmdKernel_cmdSize,
	uint64_t inscalar_writeCmdKernel_halos,
	uint64_t inscalar_writeCmdKernel_startAddress,
	const void *instream_data_in,
	size_t instream_size_data_in,
	void *outstream_data_out,
	size_t outstream_size_data_out,
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
 * \param [in] ticks_readCmdKernel0 The number of ticks for which kernel "readCmdKernel0" will run.
 * \param [in] ticks_readCmdKernel1 The number of ticks for which kernel "readCmdKernel1" will run.
 * \param [in] ticks_readCmdKernel2 The number of ticks for which kernel "readCmdKernel2" will run.
 * \param [in] ticks_writeCmdKernel The number of ticks for which kernel "writeCmdKernel" will run.
 * \param [in] inscalar_diracKernel_alpha Input scalar parameter "diracKernel.alpha".
 * \param [in] inscalar_diracKernel_beta_s Input scalar parameter "diracKernel.beta_s".
 * \param [in] inscalar_diracKernel_beta_t_b Input scalar parameter "diracKernel.beta_t_b".
 * \param [in] inscalar_diracKernel_beta_t_f Input scalar parameter "diracKernel.beta_t_f".
 * \param [in] inscalar_diracKernel_doSub Input scalar parameter "diracKernel.doSub".
 * \param [in] inscalar_diracKernel_ieo Input scalar parameter "diracKernel.ieo".
 * \param [in] inscalar_diracKernel_isign Input scalar parameter "diracKernel.isign".
 * \param [in] inscalar_readCmdKernel0_burstsPerTSlice Input scalar parameter "readCmdKernel0.burstsPerTSlice".
 * \param [in] inscalar_readCmdKernel0_cmdSize Input scalar parameter "readCmdKernel0.cmdSize".
 * \param [in] inscalar_readCmdKernel0_halos Input scalar parameter "readCmdKernel0.halos".
 * \param [in] inscalar_readCmdKernel0_startAddress Input scalar parameter "readCmdKernel0.startAddress".
 * \param [in] inscalar_readCmdKernel1_burstsPerTSlice Input scalar parameter "readCmdKernel1.burstsPerTSlice".
 * \param [in] inscalar_readCmdKernel1_cmdSize Input scalar parameter "readCmdKernel1.cmdSize".
 * \param [in] inscalar_readCmdKernel1_halos Input scalar parameter "readCmdKernel1.halos".
 * \param [in] inscalar_readCmdKernel1_startAddress Input scalar parameter "readCmdKernel1.startAddress".
 * \param [in] inscalar_readCmdKernel2_burstsPerTSlice Input scalar parameter "readCmdKernel2.burstsPerTSlice".
 * \param [in] inscalar_readCmdKernel2_cmdSize Input scalar parameter "readCmdKernel2.cmdSize".
 * \param [in] inscalar_readCmdKernel2_halos Input scalar parameter "readCmdKernel2.halos".
 * \param [in] inscalar_readCmdKernel2_startAddress Input scalar parameter "readCmdKernel2.startAddress".
 * \param [in] inscalar_writeCmdKernel_burstsPerTSlice Input scalar parameter "writeCmdKernel.burstsPerTSlice".
 * \param [in] inscalar_writeCmdKernel_cmdSize Input scalar parameter "writeCmdKernel.cmdSize".
 * \param [in] inscalar_writeCmdKernel_halos Input scalar parameter "writeCmdKernel.halos".
 * \param [in] inscalar_writeCmdKernel_startAddress Input scalar parameter "writeCmdKernel.startAddress".
 * \param [in] instream_data_in Stream "data_in".
 * \param [in] instream_size_data_in The size of the stream instream_data_in in bytes.
 * \param [out] outstream_data_out Stream "data_out".
 * \param [in] outstream_size_data_out The size of the stream outstream_data_out in bytes.
 * \param [in] routing_string A string containing comma-separated "from_name -> to_name" routing commands.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *LQCD_nonblock(
	uint64_t ticks_diracKernel,
	uint64_t ticks_readCmdKernel0,
	uint64_t ticks_readCmdKernel1,
	uint64_t ticks_readCmdKernel2,
	uint64_t ticks_writeCmdKernel,
	double inscalar_diracKernel_alpha,
	double inscalar_diracKernel_beta_s,
	double inscalar_diracKernel_beta_t_b,
	double inscalar_diracKernel_beta_t_f,
	uint64_t inscalar_diracKernel_doSub,
	uint64_t inscalar_diracKernel_ieo,
	uint64_t inscalar_diracKernel_isign,
	uint64_t inscalar_readCmdKernel0_burstsPerTSlice,
	uint64_t inscalar_readCmdKernel0_cmdSize,
	uint64_t inscalar_readCmdKernel0_halos,
	uint64_t inscalar_readCmdKernel0_startAddress,
	uint64_t inscalar_readCmdKernel1_burstsPerTSlice,
	uint64_t inscalar_readCmdKernel1_cmdSize,
	uint64_t inscalar_readCmdKernel1_halos,
	uint64_t inscalar_readCmdKernel1_startAddress,
	uint64_t inscalar_readCmdKernel2_burstsPerTSlice,
	uint64_t inscalar_readCmdKernel2_cmdSize,
	uint64_t inscalar_readCmdKernel2_halos,
	uint64_t inscalar_readCmdKernel2_startAddress,
	uint64_t inscalar_writeCmdKernel_burstsPerTSlice,
	uint64_t inscalar_writeCmdKernel_cmdSize,
	uint64_t inscalar_writeCmdKernel_halos,
	uint64_t inscalar_writeCmdKernel_startAddress,
	const void *instream_data_in,
	size_t instream_size_data_in,
	void *outstream_data_out,
	size_t outstream_size_data_out,
	const char * routing_string);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint64_t ticks_diracKernel; /**<  [in] The number of ticks for which kernel "diracKernel" will run. */
	uint64_t ticks_readCmdKernel0; /**<  [in] The number of ticks for which kernel "readCmdKernel0" will run. */
	uint64_t ticks_readCmdKernel1; /**<  [in] The number of ticks for which kernel "readCmdKernel1" will run. */
	uint64_t ticks_readCmdKernel2; /**<  [in] The number of ticks for which kernel "readCmdKernel2" will run. */
	uint64_t ticks_writeCmdKernel; /**<  [in] The number of ticks for which kernel "writeCmdKernel" will run. */
	double inscalar_diracKernel_alpha; /**<  [in] Input scalar parameter "diracKernel.alpha". */
	double inscalar_diracKernel_beta_s; /**<  [in] Input scalar parameter "diracKernel.beta_s". */
	double inscalar_diracKernel_beta_t_b; /**<  [in] Input scalar parameter "diracKernel.beta_t_b". */
	double inscalar_diracKernel_beta_t_f; /**<  [in] Input scalar parameter "diracKernel.beta_t_f". */
	uint64_t inscalar_diracKernel_doSub; /**<  [in] Input scalar parameter "diracKernel.doSub". */
	uint64_t inscalar_diracKernel_ieo; /**<  [in] Input scalar parameter "diracKernel.ieo". */
	uint64_t inscalar_diracKernel_isign; /**<  [in] Input scalar parameter "diracKernel.isign". */
	uint64_t inscalar_readCmdKernel0_burstsPerTSlice; /**<  [in] Input scalar parameter "readCmdKernel0.burstsPerTSlice". */
	uint64_t inscalar_readCmdKernel0_cmdSize; /**<  [in] Input scalar parameter "readCmdKernel0.cmdSize". */
	uint64_t inscalar_readCmdKernel0_halos; /**<  [in] Input scalar parameter "readCmdKernel0.halos". */
	uint64_t inscalar_readCmdKernel0_startAddress; /**<  [in] Input scalar parameter "readCmdKernel0.startAddress". */
	uint64_t inscalar_readCmdKernel1_burstsPerTSlice; /**<  [in] Input scalar parameter "readCmdKernel1.burstsPerTSlice". */
	uint64_t inscalar_readCmdKernel1_cmdSize; /**<  [in] Input scalar parameter "readCmdKernel1.cmdSize". */
	uint64_t inscalar_readCmdKernel1_halos; /**<  [in] Input scalar parameter "readCmdKernel1.halos". */
	uint64_t inscalar_readCmdKernel1_startAddress; /**<  [in] Input scalar parameter "readCmdKernel1.startAddress". */
	uint64_t inscalar_readCmdKernel2_burstsPerTSlice; /**<  [in] Input scalar parameter "readCmdKernel2.burstsPerTSlice". */
	uint64_t inscalar_readCmdKernel2_cmdSize; /**<  [in] Input scalar parameter "readCmdKernel2.cmdSize". */
	uint64_t inscalar_readCmdKernel2_halos; /**<  [in] Input scalar parameter "readCmdKernel2.halos". */
	uint64_t inscalar_readCmdKernel2_startAddress; /**<  [in] Input scalar parameter "readCmdKernel2.startAddress". */
	uint64_t inscalar_writeCmdKernel_burstsPerTSlice; /**<  [in] Input scalar parameter "writeCmdKernel.burstsPerTSlice". */
	uint64_t inscalar_writeCmdKernel_cmdSize; /**<  [in] Input scalar parameter "writeCmdKernel.cmdSize". */
	uint64_t inscalar_writeCmdKernel_halos; /**<  [in] Input scalar parameter "writeCmdKernel.halos". */
	uint64_t inscalar_writeCmdKernel_startAddress; /**<  [in] Input scalar parameter "writeCmdKernel.startAddress". */
	const void *instream_data_in; /**<  [in] Stream "data_in". */
	size_t instream_size_data_in; /**<  [in] The size of the stream instream_data_in in bytes. */
	void *outstream_data_out; /**<  [out] Stream "data_out". */
	size_t outstream_size_data_out; /**<  [in] The size of the stream outstream_data_out in bytes. */
	const char * routing_string; /**<  [in] A string containing comma-separated "from_name -> to_name" routing commands. */
} LQCD_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void LQCD_run(
	max_engine_t *engine,
	LQCD_actions_t *interface_actions);

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
max_run_t *LQCD_run_nonblock(
	max_engine_t *engine,
	LQCD_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void LQCD_run_group(max_group_t *group, LQCD_actions_t *interface_actions);

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
max_run_t *LQCD_run_group_nonblock(max_group_t *group, LQCD_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void LQCD_run_array(max_engarray_t *engarray, LQCD_actions_t *interface_actions[]);

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
max_run_t *LQCD_run_array_nonblock(max_engarray_t *engarray, LQCD_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* LQCD_convert(max_file_t *maxfile, LQCD_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* LQCD_init(void);

/* Error handling functions */
int LQCD_has_errors(void);
const char* LQCD_get_errors(void);
void LQCD_clear_errors(void);
/* Free statically allocated maxfile data */
void LQCD_free(void);
/* returns: -1 = error running command; 0 = no error reported */
int LQCD_simulator_start(void);
/* returns: -1 = error running command; 0 = no error reported */
int LQCD_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_LQCD_H */

