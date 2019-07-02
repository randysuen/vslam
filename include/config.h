//
// Created by randy on 19-6-14.
//

#ifndef VSLAM_CONFIG_H
#define VSLAM_CONFIG_H
#include "common_include.h"

namespace vslam
{
    class Config
    {
    private:
        static std::shared_ptr<Config> config_;
        cv::FileStorage file_;

        Config () {} // private constructor makes a singleton
    public:
        ~Config();  // close the file when deconstructing

        // set a new config file
        static void setParameterFile( const std::string& filename );

        // access the parameter values
        template< typename T >
        static T get( const std::string& key )
        {
            return T( Config::config_->file_[key] );
        }
    };
}



#endif
