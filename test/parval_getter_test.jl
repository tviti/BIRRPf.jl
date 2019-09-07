using Test
import BIRRPf


# We really just need to make sure that none of these throw an exception. Also,
# none of the values in parameters.h should ever be negative, so we can use that
# as a naive check against incorrect retvals
macro tester(varname)
    return Meta.parse("println(@test BIRRPf.get_"*string(varname)*" > 0)")
end


function main()
    @testset "parvel_getter_test" begin
        @tester "npcsm"
        @tester "nptsm"
        @tester "nptssm"
        @tester "noutm"
        @tester "ninpm"
        @tester "nrefm"
        @tester "nrsitem"
        @tester "nsectm"
    end
end


main()
