##############################################################
# Julia type for Mori dream spaces that are toric varieties
# #############################################################
abstract type ToricVarietyMDS <: MoriDreamSpace end

is_toric(::ToricVarietyMDS) = true

cox_ring_relations(::ToricVarietyMDS) = []


