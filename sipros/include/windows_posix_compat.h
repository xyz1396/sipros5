#ifndef SIPROS_WINDOWS_POSIX_COMPAT_H
#define SIPROS_WINDOWS_POSIX_COMPAT_H

#ifdef _WIN32

#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <fcntl.h>
#include <io.h>
#include <mutex>
#include <process.h>
#include <sys/stat.h>

// This adapter deliberately calls the Windows CRT's POSIX-compatible entry
// points. Keep that implementation detail from producing warnings at every
// call site while retaining warnings for the rest of the project.
#pragma warning(push)
#pragma warning(disable : 4996)

using ssize_t = std::intptr_t;

#ifndef O_CLOEXEC
#define O_CLOEXEC _O_NOINHERIT
#endif
#ifndef PROT_READ
#define PROT_READ 1
#endif
#ifndef MAP_SHARED
#define MAP_SHARED 1
#endif
#ifndef MAP_FAILED
#define MAP_FAILED reinterpret_cast<void *>(-1)
#endif

// Use the 64-bit CRT file metadata type where Unix code spells `struct stat`.
#define stat _stat64
#define off_t __int64

inline int sipros_getpid()
{
	return _getpid();
}

#define getpid sipros_getpid

inline int sipros_open(const char *path, int flags)
{
	return _open(path, flags | _O_BINARY);
}

inline int sipros_open(const char *path, int flags, int mode)
{
	return _open(path, flags | _O_BINARY, mode);
}

inline int sipros_open(const wchar_t *path, int flags)
{
	return _wopen(path, flags | _O_BINARY);
}

inline int sipros_open(const wchar_t *path, int flags, int mode)
{
	return _wopen(path, flags | _O_BINARY, mode);
}

#define open sipros_open

inline int sipros_close(int fd)
{
	return _close(fd);
}

inline int fstat(int fd, struct _stat64 *status)
{
	return _fstat64(fd, status);
}

inline int ftruncate(int fd, __int64 size)
{
	return _chsize_s(fd, static_cast<__int64>(size)) == 0 ? 0 : -1;
}

inline ssize_t pwrite(int fd, const void *buffer, size_t count, __int64 offset)
{
	// MSVC's CRT has no pwrite. Serializing the seek/write pair preserves the
	// positional-write semantics used by the parallel index builder.
	static std::mutex writeMutex;
	std::lock_guard<std::mutex> lock(writeMutex);
	if (_lseeki64(fd, offset, SEEK_SET) < 0)
		return -1;
	const unsigned int chunk = static_cast<unsigned int>(
		count > static_cast<size_t>(UINT_MAX) ? UINT_MAX : count);
	return static_cast<ssize_t>(_write(fd, buffer, chunk));
}

inline void *mmap(void *, size_t size, int, int, int fd, __int64)
{
	const intptr_t osHandle = _get_osfhandle(fd);
	if (osHandle == -1)
	{
		errno = EBADF;
		return MAP_FAILED;
	}
	HANDLE mapping = CreateFileMappingW(
		reinterpret_cast<HANDLE>(osHandle), nullptr, PAGE_READONLY, 0, 0, nullptr);
	if (mapping == nullptr)
	{
		errno = EIO;
		return MAP_FAILED;
	}
	void *view = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, size);
	CloseHandle(mapping);
	if (view == nullptr)
	{
		errno = EIO;
		return MAP_FAILED;
	}
	return view;
}

inline int munmap(void *address, size_t)
{
	if (UnmapViewOfFile(address) == 0)
	{
		errno = EIO;
		return -1;
	}
	return 0;
}

#pragma warning(pop)

#endif // _WIN32

#endif // SIPROS_WINDOWS_POSIX_COMPAT_H
