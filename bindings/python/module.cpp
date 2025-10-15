#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <optional>
#include <nlohmann/json.hpp>
#include "../../core/registry.h"
#include "../../include/Matrix.h"

extern "C" void force_link_adapters();

namespace py = pybind11;
using namespace core;
using my_matrix::Matrix;

static py::object json_to_py(const nlohmann::json& j) {
    auto json = py::module_::import("json");
    return json.attr("loads")(py::str(j.dump()));
}


static Matrix matrix_from_numpy(const py::array_t<double, py::array::c_style | py::array::forcecast>& arr) {
    if (arr.ndim() != 2) throw std::runtime_error("matrix must be 2D");
    const unsigned int rows = static_cast<unsigned int>(arr.shape(0));
    const unsigned int cols = static_cast<unsigned int>(arr.shape(1));
    Matrix M(rows, cols);  

    auto r = arr.unchecked<2>();
    for (unsigned int i = 0; i < rows; ++i) {
        for (unsigned int j = 0; j < cols; ++j) {
            M[i][j] = r(i, j);
        }
    }
    return M;
}

PYBIND11_MODULE(uav_py, m) {
    m.doc() = "Python bindings for UAV core";
    force_link_adapters();
    py::class_<Result>(m, "Result")
        .def_readonly("path", &Result::path)
        .def_readonly("cost", &Result::cost)
        .def_readonly("algorithm", &Result::algorithm)
        .def_readonly("problem", &Result::problem)
        .def_property_readonly(
            "meta",
            [](const Result& r) {
                auto json = py::module_::import("json");
                return json.attr("loads")(py::str(r.meta.dump()));
            }
        );
    m.def(
        "list_algorithms",
        [](std::optional<std::string> problem) {
            auto all = core::Registry::instance().list(problem);
            py::list out;  
            for (const auto& a : all) {
                py::dict d;
                d["key"] = a.key;
                d["problem"] = a.problem;
                d["schema"] = json_to_py(a.schema);
                out.append(std::move(d));
            }
            return out;
        },
        py::arg("problem") = py::none()
    );

    m.def(
        "run",
        [](const std::string& algo,
            const py::array_t<double, py::array::c_style | py::array::forcecast>& matrix,
            py::object params_obj /* dict | None */) {
                Matrix M = matrix_from_numpy(matrix);
                nlohmann::json params = nlohmann::json::object();
                if (!params_obj.is_none()) {
                    auto json = py::module_::import("json");
                    std::string dumped = json.attr("dumps")(params_obj).cast<std::string>();
                    params = nlohmann::json::parse(dumped);
                }

                const auto& info = Registry::instance().get(algo);
                return info.run(M, params);
        },
        py::arg("algorithm"),
        py::arg("matrix"),
        py::arg("params") = py::none()
    );
}
