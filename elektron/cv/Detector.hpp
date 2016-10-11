#ifndef ELEKTRON_DETECTOR_H_
#define ELEKTRON_DETECTOR_H_

namespace elektron
{

enum DetectorType = {canny, st}

template <typename T_, std::size_t R_, std::size_t C_>
class Detector
{
public:
	Detector();
	Detector(DetectorType d) : det_tp_(d) { };

	void detect(const Image<T_, R_, C_>& im);
	void set_type(DetectorType d) : det_tp_(d) { };
	
private:
	DetectorType det_tp_;

	void shi_tomasi();
};

template <typename T_, std::size_t R_, std::size_t C_>
void detect(const Image<T_, R_, C_>& im)
{
	switch (det_tp_)
	{
		case canny:
			break;
		case st:
			break;
		default:
			break;
	}
}

template <typename T_, std::size_t R_, std::size_t C_>
void shi_tomasi(const Image<T_, R_, C_>& im)
{

}

}
#endif
