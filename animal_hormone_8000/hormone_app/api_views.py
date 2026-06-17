from rest_framework import viewsets, permissions, status
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework.decorators import action
from django.contrib.auth.models import User
from django.contrib.auth import authenticate, login
from rest_framework_simplejwt.tokens import RefreshToken
from .models import Hormone, Literature, AnalysisHistory, AnalysisTask
from .serializers import (
    HormoneSerializer, LiteratureSerializer, 
    AnalysisHistorySerializer, UserSerializer
)
import json
import logging

logger = logging.getLogger(__name__)


class HormoneViewSet(viewsets.ModelViewSet):
    """API端点，用于激素受体数据的查看和编辑"""
    queryset = Hormone.objects.all()
    serializer_class = HormoneSerializer
    permission_classes = [permissions.IsAuthenticatedOrReadOnly]


class LiteratureViewSet(viewsets.ModelViewSet):
    """API端点，用于文献数据的查看和编辑"""
    queryset = Literature.objects.all()
    serializer_class = LiteratureSerializer
    permission_classes = [permissions.IsAuthenticatedOrReadOnly]
    
    @action(detail=True, methods=['get'])
    def download(self, request, pk=None):
        """下载文献文件的自定义动作"""
        literature = self.get_object()
        # 文件下载逻辑...


class AnalysisHistoryViewSet(viewsets.ModelViewSet):
    """API端点，用于分析历史的查看和创建"""
    queryset = AnalysisHistory.objects.all()
    serializer_class = AnalysisHistorySerializer
    permission_classes = [permissions.IsAuthenticated]
    
    def get_queryset(self):
        """限制用户只能查看自己的分析历史"""
        return self.queryset.filter(user=self.request.user)


class GeneAnalysisAPIView(APIView):
    """基因分析API端点"""
    permission_classes = [permissions.IsAuthenticated]
    
    def post(self, request):
        """处理基因序列分析请求"""
        sequence = request.data.get('sequence')
        analysis_type = request.data.get('analysis_type', 'basic')
        
        if not sequence:
            return Response(
                {'error': '序列不能为空'}, 
                status=status.HTTP_400_BAD_REQUEST
            )
        
        try:
            # 执行分析（实际项目中可能会调用分析函数或异步任务）
            results = perform_analysis(sequence, analysis_type)
            
            # 保存分析历史
            history = AnalysisHistory.objects.create(
                user=request.user,
                sequence=sequence[:10000],  # 限制序列长度
                results=results
            )
            
            return Response(
                {'results': results, 'history_id': history.id},
                status=status.HTTP_200_OK
            )
            
        except Exception as e:
            logger.error(f"分析失败: {str(e)}")
            return Response(
                {'error': '分析过程出错，请重试'},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR
            )


class UserRegistrationAPIView(APIView):
    """用户注册API"""
    permission_classes = [permissions.AllowAny]
    
    def post(self, request):
        serializer = UserSerializer(data=request.data)
        if serializer.is_valid():
            user = serializer.save()
            refresh = RefreshToken.for_user(user)
            return Response({
                'refresh': str(refresh),
                'access': str(refresh.access_token),
                'user': serializer.data
            }, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class UserLoginAPIView(APIView):
    """用户登录API"""
    permission_classes = [permissions.AllowAny]
    
    def post(self, request):
        username = request.data.get('username')
        password = request.data.get('password')
        
        user = authenticate(username=username, password=password)
        if user is not None:
            login(request, user)
            refresh = RefreshToken.for_user(user)
            return Response({
                'refresh': str(refresh),
                'access': str(refresh.access_token),
                'user': UserSerializer(user).data
            })
        return Response({'error': '无效的凭据'}, status=status.HTTP_401_UNAUTHORIZED)


class CurrentUserView(APIView):
    """获取当前用户信息"""
    permission_classes = [permissions.IsAuthenticated]
    
    def get(self, request):
        serializer = UserSerializer(request.user)
        return Response(serializer.data)


class TaskListView(APIView):
    """异步任务列表API"""
    permission_classes = [permissions.IsAuthenticated]
    
    def get(self, request):
        tasks = AnalysisTask.objects.filter(user=request.user)
        # 返回任务列表...


class TaskDetailView(APIView):
    """异步任务详情API"""
    permission_classes = [permissions.IsAuthenticated]
    
    def get(self, request, task_id):
        try:
            task = AnalysisTask.objects.get(task_id=task_id, user=request.user)
            # 返回任务详情...
        except AnalysisTask.DoesNotExist:
            return Response({'error': '任务不存在'}, status=404)    