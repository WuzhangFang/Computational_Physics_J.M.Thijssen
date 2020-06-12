#ifndef ASSERT_
#define ASSERT_

#ifdef NO_ARG_CHECK
#define Assert(condition, message)
#else /* NO_ARG_CHECK */
#include <iostream>
#define Assert(condition, message)\
{\
  if(!(condition)) std::cerr << (message) << std::endl;\
}
#endif /* NO_ARG_CHECK */

#ifdef T_DEBUG
#define T_LOG(x) x
#else
#define T_LOG(x)
#endif

#endif /* ASSERT_ */
