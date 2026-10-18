/*************************************************************************
    > File Name: Msg.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Fri 31 Jul 2015 03:20:14 PM CDT
 ************************************************************************/

#include "utility/src/msg.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

DREAMPLACE_BEGIN_NAMESPACE

int dreamplacePrint(MessageType m, const char* format, ...) {
  va_list args;
  va_start(args, format);
  int ret = dreamplaceVPrintStream(m, stdout, format, args);
  va_end(args);

  return ret;
}

int dreamplacePrintStream(MessageType m, FILE* stream, const char* format,
                          ...) {
  va_list args;
  va_start(args, format);
  int ret = dreamplaceVPrintStream(m, stream, format, args);
  va_end(args);

  return ret;
}

int dreamplaceVPrintStream(MessageType m, FILE* stream, const char* format,
                           va_list args) {
#ifdef DREAMPLACE_TORCH_FREE
  // An embedded library must not flood the host's normal output. Keep
  // diagnostics; the native API returns its placement metrics to the caller.
  if (m == kINFO || m == kDEBUG || m == kNONE) return 0;
#endif
  // print prefix
  char prefix[16];
  dreamplaceSPrintPrefix(m, prefix);
  fprintf(stream, "%s", prefix);

  // print message
  int ret = vfprintf(stream, format, args);

  return ret;
}

int dreamplaceSPrint(MessageType m, char* buf, const char* format, ...) {
  va_list args;
  va_start(args, format);
  int ret = dreamplaceVSPrint(m, buf, format, args);
  va_end(args);

  return ret;
}

int dreamplaceVSPrint(MessageType m, char* buf, const char* format,
                      va_list args) {
  // print prefix
  char prefix[16];
  dreamplaceSPrintPrefix(m, prefix);
  const size_t prefix_length = strlen(prefix);
  memcpy(buf, prefix, prefix_length + 1);

  // This legacy API requires a caller-provided buffer large enough for the
  // prefix and message. Measure using a copy so args remains available.
  va_list measure_args;
  va_copy(measure_args, args);
  const int length = vsnprintf(nullptr, 0, format, measure_args);
  va_end(measure_args);
  if (length < 0) return length;
  return vsnprintf(buf + prefix_length, static_cast<size_t>(length) + 1,
                   format, args);
}

int dreamplaceSPrintPrefix(MessageType m, char* prefix) {
  switch (m) {
    case kNONE:
      prefix[0] = prefix[1] = '\0';
      return 1;
    case kINFO:
      memcpy(prefix, "[INFO   ] ", 11);
      return 10;
    case kWARN:
      memcpy(prefix, "[WARNING] ", 11);
      return 10;
    case kERROR:
      memcpy(prefix, "[ERROR  ] ", 11);
      return 10;
    case kDEBUG:
      memcpy(prefix, "[DEBUG  ] ", 11);
      return 10;
    case kASSERT:
      memcpy(prefix, "[ASSERT ] ", 11);
      return 10;
    default:
      dreamplaceAssertMsg(0, "unknown message type");
  }
  return 0;
}

void dreamplacePrintAssertMsg(const char* expr, const char* fileName,
                              unsigned lineNum, const char* funcName,
                              const char* format, ...) {
  // construct message
  char buf[1024];
  va_list args;
  va_start(args, format);
  vsnprintf(buf, sizeof(buf), format, args);
  va_end(args);

  // print message
  dreamplacePrintStream(kASSERT, stderr,
                        "%s:%u: %s: Assertion `%s' failed: %s\n", fileName,
                        lineNum, funcName, expr, buf);
}

void dreamplacePrintAssertMsg(const char* expr, const char* fileName,
                              unsigned lineNum, const char* funcName) {
  // print message
  dreamplacePrintStream(kASSERT, stderr, "%s:%u: %s: Assertion `%s' failed\n",
                        fileName, lineNum, funcName, expr);
}

DREAMPLACE_END_NAMESPACE
