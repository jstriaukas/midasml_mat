function out = standardize(in)
m_  = mean(in);
s_  = std(in);
out = (in - m_)./s_;
end
