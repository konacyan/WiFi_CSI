function  perm( list,s,e )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
global matrix_perm perm_idx;

    if(s == e)
%         disp(num2str(list));
        matrix_perm(perm_idx,:) = list;
        perm_idx = perm_idx + 1;

    else
        for (i = s:e)
            t = list(s);list(s) = list(i);list(i) = t;
            perm(list,s+1,e);
            t = list(s);list(s) = list(i);list(i) = t;
        end
    end
end

