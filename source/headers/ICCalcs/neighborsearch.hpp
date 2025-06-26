#pragma once

#include <vector>
#include <array>
#include <utility>
#include <cmath>
#include <memory>
#include <iostream>
#include <algorithm>
#include <omp.h>

#include "external/nanoflann/nanoflann.hpp"

class NeighborSearch {
private:
    struct PointCloud {
        std::vector<std::array<double, 3>> points;

        inline size_t kdtree_get_point_count() const {
            return points.size();
        }

        inline double kdtree_get_pt(const size_t idx, const size_t axis) const {
            return points[idx][axis];
        }

        template <class BBOX>
        bool kdtree_get_bbox(BBOX&) const {
            return false;
        }
    };

    using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, PointCloud>,
        PointCloud,
        3, uint32_t
    >;

    PointCloud cloud;
    std::unique_ptr<KDTree> index;

public:
    NeighborSearch() = default;
    ~NeighborSearch() = default;

    // Build KDTree from the set of points
    void build(const std::vector<std::array<double, 3>>& points) {
        cloud.points = points;
        index.reset(new KDTree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
        index->buildIndex();
    }

    // We replace std::pair with this one for performance and also avoid sqrt in query function
    struct Contact {
        int ligand_idx;
        int protein_idx;
        double distance_squared;
    };

    // Query neighbors within the interval distance [min, max] using RadiusResultSet and OpenMP
    std::vector<Contact> query(const std::vector<std::array<double, 3>>& queryPoints,
                               double dist_max, double dist_min = 0.0) const {
        std::vector<Contact> result;
        if (!index) 
            return result;

        const double max_radius_sq = dist_max * dist_max;
        const double min_radius_sq = dist_min * dist_min;

        size_t guess = std::min<size_t>(cloud.points.size(), 8);
        result.reserve(queryPoints.size() * guess);

        #pragma omp parallel
        {
            std::vector<Contact> local_result;
            std::vector<nanoflann::ResultItem<uint32_t, double>> local_matches;

            #pragma omp for nowait
            for (int i = 0; i < static_cast<int>(queryPoints.size()); ++i) {
                local_matches.clear();
                nanoflann::RadiusResultSet<double, uint32_t> resultSet(max_radius_sq, local_matches);
                const double* query_point = queryPoints[i].data();
                nanoflann::SearchParameters params;
                index->findNeighbors(resultSet, query_point, params);

                for (const auto& match : local_matches) {
                    if (match.second >= min_radius_sq) {
                        local_result.push_back({ i, static_cast<int>(match.first), match.second });
                    }
                }
            }

            #pragma omp critical
            result.insert(result.end(), local_result.begin(), local_result.end());
        }

        return std::move(result);
    }
};
