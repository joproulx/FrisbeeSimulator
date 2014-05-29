using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Windows.Data;
using System.Windows.Media;

namespace FlightSimulator.Converters
{
    [ValueConversion(typeof(double), typeof(double))]
    public class DivisionConverter : IValueConverter
    {
        public object Convert(object value, Type targetType,
            object parameter, CultureInfo culture)
        {
            double dividend = (double)value;
            double divisor = double.Parse((string)parameter);

            return (double)dividend / divisor;
        }

        public object ConvertBack(object value, Type targetType,
            object parameter, CultureInfo culture)
        {
            return null;
        }
    }
}
