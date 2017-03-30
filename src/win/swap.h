
#pragma once

#ifdef _DSO_ON_WIN

namespace dso
{

	template <typename t>
	void swap(t& x, t& y)
	{
		t temp = x;
		x = y;
		y = temp;
		return;
	}

}

#endif //_DSO_ON_WIN