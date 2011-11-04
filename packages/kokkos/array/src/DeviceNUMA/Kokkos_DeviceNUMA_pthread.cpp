/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2011 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

/*--------------------------------------------------------------------------*/
/* Kokkos interfaces */

#include <Kokkos_DeviceNUMA.hpp>

/*--------------------------------------------------------------------------*/
/* Standard 'C' Linux libraries */
#include <pthread.h>
#include <sched.h>
#include <errno.h>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// Driver for each created pthread

namespace {

void * device_numa_pthread_driver( void * arg )
{
  ((DeviceNUMAThread *) arg)->driver();
  return NULL ;
}

pthread_mutex_t device_numa_pthread_mutex = PTHREAD_MUTEX_INITIALIZER ;

}

//----------------------------------------------------------------------------
// Spawn this thread

bool device_numa_thread_spawn( DeviceNUMAThread * thread )
{
  bool result = false ;

  pthread_attr_t attr ;
  
  if ( 0 == pthread_attr_init( & attr ) ||
       0 == pthread_attr_setscope(       & attr, PTHREAD_SCOPE_SYSTEM ) ||
       0 == pthread_attr_setdetachstate( & attr, PTHREAD_CREATE_DETACHED ) ) {

    pthread_t pt ;

    result =
      0 == pthread_create( & pt, & attr, device_numa_pthread_driver, thread );
  }

  pthread_attr_destroy( & attr );

  return result ;
}

//----------------------------------------------------------------------------
// Mutually exclusive locking and unlocking

void device_numa_thread_lock()
{ pthread_mutex_lock( & device_numa_pthread_mutex ); }

void device_numa_thread_unlock()
{ pthread_mutex_unlock( & device_numa_pthread_mutex ); }

//----------------------------------------------------------------------------
// Performance critical function: thread waits while value == *state

void DeviceNUMAThread::wait( const DeviceNUMAThread::State flag )
{
  const long value = flag ;
  while ( value == m_state ) {
    sched_yield();
  }
}

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

